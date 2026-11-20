"""Normalize existing covariance tables; never impute cross-programme zeros."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent


def assemble():
    catalogue = pd.read_csv(ROOT / "effect_catalogue.csv")
    assert catalogue.effect_id.is_unique
    variance = catalogue.set_index("effect_id").sampling_variance
    rows, blocks, hashes = [], [], {}
    for source, effects in catalogue.groupby("covariance_table", sort=False):
        path = ROOT / source
        hashes[source] = hashlib.sha256(path.read_bytes()).hexdigest()
        table = pd.read_csv(path)
        ids = effects.effect_id.tolist()
        assert effects.programme_id.nunique() == 1
        programme = effects.programme_id.iloc[0]
        if "effect_id_1" in table:
            table = table.rename(columns={"covariance": "sampling_covariance"})
        elif "effect_id_i" in table:
            table = table.rename(columns={"effect_id_i": "effect_id_1", "effect_id_j": "effect_id_2", "covariance": "sampling_covariance"})
        else:
            # The two wide tables use endpoint labels, not global effect IDs.
            assert effects.endpoint.is_unique
            lookup = dict(zip(effects.endpoint, effects.effect_id))
            index = table.columns[0]
            assert set(table[index]) == set(lookup)
            assert set(table.columns[1:]) == set(lookup)
            table = table.melt(id_vars=index, var_name="second", value_name="sampling_covariance")
            table["effect_id_1"] = table[index].map(lookup)
            table["effect_id_2"] = table.second.map(lookup)
        table = table.loc[table.effect_id_1.isin(ids) & table.effect_id_2.isin(ids)].copy()
        assert not table.duplicated(["effect_id_1", "effect_id_2"]).any()
        assert len(table) == len(ids) ** 2, source
        matrix = table.pivot(index="effect_id_1", columns="effect_id_2", values="sampling_covariance").reindex(index=ids, columns=ids).to_numpy(float)
        assert np.isfinite(matrix).all()
        assert np.allclose(matrix, matrix.T, rtol=1e-10, atol=1e-12), source
        assert np.allclose(np.diag(matrix), variance.loc[ids], rtol=1e-10, atol=1e-12), source
        eigmin = float(np.linalg.eigvalsh((matrix + matrix.T) / 2).min())
        assert eigmin >= -1e-9 * max(1.0, float(np.abs(matrix).max())), (source, eigmin)
        table["programme_id"] = programme
        table["source_covariance_table"] = source
        table["status"] = "source_covariance_not_final_adjudication"
        rows.append(table[["effect_id_1", "effect_id_2", "sampling_covariance", "programme_id", "source_covariance_table", "status"]])
        blocks.append(dict(programme_id=programme, effects=len(ids), covariance_entries=len(table), minimum_eigenvalue=eigmin))
    combined = pd.concat(rows, ignore_index=True)
    assert set(combined.effect_id_1) == set(catalogue.effect_id)
    assert not combined.duplicated(["effect_id_1", "effect_id_2"]).any()
    combined.to_csv(ROOT / "catalogue_covariance.csv", index=False)
    hashes["effect_catalogue.csv"] = hashlib.sha256((ROOT / "effect_catalogue.csv").read_bytes()).hexdigest()
    summary = dict(effect_rows=len(catalogue), source_blocks=len(blocks), covariance_entries=len(combined),
                   unrepresented_ordered_pairs=len(catalogue)**2-len(combined),
                   missing_cross_block_pairs="unknown_not_zero; independence_not_adjudicated",
                   blocks=blocks, input_sha256=hashes,
                   limitation="Source design assumptions retained. Numerical PSD does not validate sampling independence. No pooled estimate or eligibility decision added.")
    (ROOT / "catalogue_covariance_audit.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({k:v for k,v in summary.items() if k not in ["blocks", "input_sha256"]}, indent=2))


if __name__ == "__main__":
    assemble()
