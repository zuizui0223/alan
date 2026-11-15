"""Offline gate: nine original ALAN plant reproductive publications, no fabricated pooling."""
import json
from pathlib import Path
from scripts.analysis.alan_A_stage_estimand_evidence_guard_20261115 import validate, ROOT, MAP, REG
def _data():
    return (json.loads(MAP.read_text(encoding="utf-8")),
            json.loads(REG.read_text(encoding="utf-8")))
def test_nine_papers_no_new_abstract_count():
    a,b=_data(); r=validate(a,b)
    assert r["direct_manipulated_papers"]==9
    assert r["additional_unverified_conference_leads"]==1
    assert r["confirmed_ALAN_x_independent_dark_V_meta_effects"]==0
def test_unknown_lu_open_seed_never_becomes_zero():
    a,b=_data()
    lu=next(x for x in a["papers"] if x["id"]=="lu2026jpe")
    lu["outcomes"][1]["value"]=0.0
    import pytest
    with pytest.raises(AssertionError):
        validate(a,b)
def test_abstract_not_made_an_experiment():
    a,b=_data()
    a["additional_grey_leads"][0]["status"]="peer_reviewed"
    import pytest
    with pytest.raises(AssertionError):
        validate(a,b)
def test_signs_do_not_become_joint_meta_effects():
    a,b=_data()
    a["papers"][0]["outcomes"][0]["ready"]=True
    import pytest
    with pytest.raises(AssertionError):
        validate(a,b)
