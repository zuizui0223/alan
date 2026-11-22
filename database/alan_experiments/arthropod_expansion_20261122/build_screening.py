"""Assemble the dated rapid screen, identity flags and RIS without changing history.

Python standard library only. Bibliographic details are a non-abstract projection
of the hashed original index response. Screening inputs are human-readable AI
decisions, not a second independent human screen or final eligibility decisions.
"""
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
import unicodedata

ROOT = Path(__file__).resolve().parent
BASE = ROOT.parent / "search"


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write(name, rows, fields=None):
    with (ROOT / name).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def title_key(value):
    value = unicodedata.normalize("NFKC", value).casefold()
    return re.sub(r"[^\w]", "", value)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    records = read(ROOT / "records.csv")
    unmatched = read(ROOT / "unmatched_records.csv")
    ids = {r["record_id"] for r in unmatched}
    assert len(ids) == len(unmatched)
    assert ids == {r["record_id"] for r in records if r["match_status"] == "unmatched_candidate_report"}
    source_paths = sorted((ROOT / "screening_sources").glob("*.csv"))
    decisions = {}
    for path in source_paths:
        for row in read(path):
            ident = row["record_id"]
            assert ident not in decisions, f"Duplicate screen: {ident}"
            decisions[ident] = {**row, "screening_input": str(path.relative_to(ROOT))}
    assert set(decisions) == ids, "Missing or extra screening decision"
    allowed = {"fulltext_needed", "separate_observational", "separate_indoor",
               "separate_spectrum_no_dark", "reference_mining", "out_of_scope", "identity_review"}
    screened = []
    for row in unmatched:
        decision = decisions[row["record_id"]]
        assert decision["candidate_report_key"] == row["candidate_report_key"]
        assert decision["decision"] in allowed and decision["reason"].strip()
        screened.append({**row, **decision,
                         "screening_status": "rapid_screen_complete_not_final_eligibility"})
    write("screened_unmatched_records.csv", screened)
    fulltext = [r for r in screened if r["decision"] == "fulltext_needed"]
    write("fulltext_candidates.csv", fulltext, list(screened[0]))

    # Title matches flag possible identity; they never rewrite report keys.
    existing = read(BASE / "combined_search_records.csv")
    title_index = defaultdict(set)
    for row in existing:
        if row["title"]:
            title_index[title_key(row["title"])].add(row["report_key"])
    flags, seen = [], set()

    def flag(ident, target, scope, basis, evidence):
        key = (ident, target)
        if key in seen:
            return
        seen.add(key)
        flags.append(dict(record_id=ident, target_key=target, comparison_scope=scope,
                          evidence_basis=basis, status="candidate_not_merged", evidence=evidence))

    for row in unmatched:
        for target in sorted(title_index.get(title_key(row["title"]), set())):
            flag(row["record_id"], target, "baseline", "normalized_title",
                 "Same normalized title; verify authors, journal metadata and version relation.")
    within = defaultdict(list)
    for row in unmatched:
        within[title_key(row["title"])].append(row)
    for group in within.values():
        for index, row in enumerate(group):
            for other in group[index + 1:]:
                flag(row["record_id"], other["candidate_report_key"], "new_snapshot", "normalized_title",
                     "Same normalized title; distinguish duplicate index entries from preprint/final versions.")
    manual = [
        ("PPR:PPR926215", "doi:10.1242/jeb.250443", "baseline", "Aquatic shredder preprint/final candidate identified in abstract screen."),
        ("PPR:PPR913088", "doi:10.1371/journal.pone.0329449", "baseline", "Gammarid preprint/final candidate identified in abstract screen."),
        ("PPR:PPR1196034", "AGR:IND609294764", "new_snapshot", "Keystoneness preprint/final candidate; verify version relation."),
        ("PPR:PPR381597", "AGR:IND607762636", "new_snapshot", "Matthew-effect report/preprint candidate; verify title and authors.")
    ]
    baseline_keys = {r["report_key"] for r in read(BASE / "unified_review_inventory.csv")}
    new_keys = {r["candidate_report_key"] for r in unmatched}
    for ident, target, scope, evidence in manual:
        assert ident in ids
        assert target in (baseline_keys if scope == "baseline" else new_keys)
        flag(ident, target, scope, "screening_version_candidate", evidence)
    write("identity_candidates.csv", flags)

    details = {r["record_id"]: r for r in read(ROOT / "bibliographic_details.csv")}
    assert set(details) == {r["record_id"] for r in records}
    ris = []
    for row in unmatched:
        detail = details[row["record_id"]]
        source = row["record_id"].split(":", 1)[0]
        ris_type = {"PPR": "UNPB", "ETH": "THES", "PAT": "PAT"}.get(source, "JOUR")
        item = [("TY", ris_type),
                ("ID", row["record_id"]), ("TI", row["title"])]
        item += [("AU", name) for name in json.loads(detail["authors_json"])]
        item += [("PY", row["year"]), ("JO", row["journal"]), ("VL", detail["volume"]),
                 ("IS", detail["issue"]), ("SP", detail["pages"]),
                 ("DO", row["doi"]), ("UR", row["source_url"]),
                 ("DB", "Europe PMC; supplementary search 2026-11-22"),
                 ("N1", "Rapid screen: " + decisions[row["record_id"]]["decision"] + "; identity and eligibility not final.")]
        for tag, value in item:
            if value:
                ris.append(f"{tag}  - {' '.join(str(value).split())}")
        ris.extend(["ER  -", ""])
    (ROOT / "unmatched_records.ris").write_text("\n".join(ris), encoding="utf-8")
    assert sum(line.startswith("TY  -") for line in ris) == len(unmatched)
    counts = dict(sorted(Counter(r["decision"] for r in screened).items()))
    inputs = [ROOT / "records.csv", ROOT / "unmatched_records.csv", ROOT / "bibliographic_details.csv",
              BASE / "combined_search_records.csv", BASE / "unified_review_inventory.csv", *source_paths]
    outputs = ["screened_unmatched_records.csv", "fulltext_candidates.csv", "identity_candidates.csv", "unmatched_records.ris"]
    summary = dict(date_jst="2026-11-22", indexed_records=len(records), unmatched_records=len(unmatched),
                   screened_records=len(screened), screening_decisions=counts,
                   identity_candidate_pairs=len(flags), identity_pair_scopes=dict(Counter(r["comparison_scope"] for r in flags)),
                   duplicate_identity_candidates_merged=0, final_eligibility=False,
                   independent_experiment_count=None,
                   screening_process="AI-assisted rapid title/abstract screen, split into two non-overlapping batches; no independent human double screen.",
                   interpretation="Categories route records for follow-up; they are not final exclusions, distinct reports, or experimental programmes.",
                   input_sha256={str(p.relative_to(ROOT.parent)): sha(p) for p in inputs},
                   output_sha256={name: sha(ROOT / name) for name in outputs})
    (ROOT / "screening_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: summary[k] for k in ["indexed_records", "unmatched_records", "screening_decisions", "identity_candidate_pairs", "identity_pair_scopes"]}, indent=2))


if __name__ == "__main__":
    main()
