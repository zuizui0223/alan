#!/usr/bin/env python3
"""Original 92 published diel pollinator-observation studies: public DOI and OA gate.

This is paper-level *bibliographic acquisition*, NOT field pollinator identification.
All true field reproduction source-pairs come from existing SHA-verified 301-case
pre-exposure original paper/site/year/season source (without exposing these details).
Author 1652 treatment records are Git-blob verified. Do not infer species or
pollinator guilds from paper titles, abstracts, keywords or OA metadata.
DOI URL, public paper title, open-access candidate landing URL and source study
reference can be published; source plant/site/year identifiers cannot.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

import pandas as pd

ORIGINAL_BLOB = "39a82d966afcdde670603add1c881ace8791ffe9"
PREFLIGHT_SHA = "21b1fef1beb9e6238cfa47101b991815d4c0f33787b9737f834426dc3bab15af"
DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$", re.I)
URL_RE = re.compile(r"^https?://", re.I)


def require(x, msg):
    if not bool(x):
        raise ValueError("ALAN_ORIGINAL_PAPER_OA_GATE: " + msg)


def file_sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def gitblob(path):
    data = Path(path).read_bytes()
    return hashlib.sha1(b"blob " + str(len(data)).encode() + b"\x00" + data).hexdigest()


def normalize_identifier(raw):
    """Preserve non-DOI and missing-source statuses; never invent DOI."""
    q = str(raw or "").strip()
    if q in {"", "NA"}:
        return "missing", None
    q = re.sub(r"^doi\s*:\s*", "", q, flags=re.I)
    q = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", q, flags=re.I)
    q = urllib.parse.unquote(q).strip()
    if DOI_RE.fullmatch(q):
        return "doi_syntax", q
    if URL_RE.match(q):
        return "non_doi_url", q
    return "invalid_literal", q


def source_queue(author, preflight):
    require(gitblob(author) == ORIGINAL_BLOB, "Author study table Git blob changed")
    require(file_sha256(preflight) == PREFLIGHT_SHA,
            "Original 301 source compressed bytes changed")
    x = pd.read_csv(author, dtype=str, keep_default_na=False)
    require(len(x) == 1652 and x.study_ID.nunique() == 140,
            "Original published source 1652/140 changed")
    y = pd.read_csv(preflight, dtype=str)
    require((len(y), y.study_key.nunique(), y.species.nunique()) ==
            (301, 118, 123), "True experiment site/year cohort is not 301/118/123")
    require(y.source_pair_id.is_unique, "Source-pair key duplicated")
    paper = {}
    for key, rows in x.groupby("study_ID", observed=True):
        values = set(rows.Pollinator_observations_included)
        require(values <= {"y", "n", "NA"}, "Unrecognized study-visitor flag")
        status = ("included_some_observations" if "y" in values else
                  "not_included" if values == {"n"} else "unknown_or_mixed")
        literals = sorted(set(rows.DOI) - {"", "NA"})
        kind, identifier = normalize_identifier(literals[0] if len(literals) == 1 else "")
        if len(literals) > 1:
            kind, identifier = "ambiguous_identifier", None
        pubyear = sorted(set(rows.publication_year) - {"", "NA"})
        paper[str(key)] = {
            "visitor_observation_flag": status,
            "original_DOI_literal_status": kind,
            "original_bibliographic_identifier": identifier,
            "original_study_reference": str(key),  # published paper ID only
            "author_publication_year": (pubyear[0] if len(pubyear) == 1 else None),
            "bib_literal_count": len(literals),
        }
    counts = Counter(v["visitor_observation_flag"] for v in paper.values())
    require(counts == {
        "included_some_observations": 111, "not_included": 27, "unknown_or_mixed": 2,
    }, "Original 140 paper visitor flags changed")
    ids = set(y.study_key.astype(str))
    require(ids <= set(paper), "Original real reproductive study key unmatched")
    keys = sorted(k for k in ids if paper[k]["visitor_observation_flag"] ==
                  "included_some_observations")
    require(len(keys) == 92 and int(y.study_key.isin(keys).sum()) == 247 and
            y.loc[y.study_key.isin(keys), "species"].nunique() == 110,
            "Original 92/247/110 actual published reproduction coverage altered")
    pubs = y.study_key.value_counts().to_dict()
    result = [{
        **paper[k],
        "actual_authentic_field_daynight_reproductive_comparisons": int(pubs[k]),
        "actual_local_field_visitor_taxon_verified": False,
        "same_experiment_plant_site_year_visitor_evidence_verified": False,
    } for k in keys]
    kinds = Counter(v["original_DOI_literal_status"] for v in result)
    require(kinds == {"doi_syntax": 87, "invalid_literal": 3,
                      "non_doi_url": 1, "missing": 1},
            "True original DOI versus alternative-identifier denominator changed")
    return result, kinds


def request_json(url, timeout=12):
    # Bounded public metadata requests only: never fetch paywalled PDF or
    # copyrighted article body; don't leak original plant/experiment sites.
    req = urllib.request.Request(url, headers={
        "User-Agent": "ALAN-metadata-audit/0.2 (+https://github.com/zuizui0223/alan)",
        "Accept": "application/json",
    })
    for attempt in range(2):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as f:
                return json.load(f), "ok"
        except urllib.error.HTTPError as e:
            if e.code == 429:
                if attempt == 0:
                    time.sleep(1.5)
                    continue
                return None, "rate_limited"
            if e.code == 404:
                return None, "not_found"
            return None, "http_" + str(e.code)
        except (urllib.error.URLError, TimeoutError, OSError,
                ValueError, json.JSONDecodeError):
            if attempt == 0:
                time.sleep(.5)
                continue
            return None, "unavailable"
    return None, "unavailable"


def year_of_crossref(message):
    for k in ("published", "published-online", "issued", "published-print"):
        v = message.get(k, {})
        parts = v.get("date-parts", []) if isinstance(v, dict) else []
        if parts and parts[0] and isinstance(parts[0][0], int):
            return parts[0][0]
    return None


def metadata_for(doi, fetch=request_json):
    q = urllib.parse.quote(doi, safe="/")
    cross, cs = fetch("https://api.crossref.org/works/" + q)
    cr = cross.get("message", {}) if isinstance(cross, dict) else {}
    if not isinstance(cr, dict):
        cr = {}
    rdoi = cr.get("DOI")
    if rdoi is not None and rdoi.strip().lower() != doi.lower():
        cs = "wrong_doi"
        cr = {}
    titles = cr.get("title") or []
    title = titles[0].strip() if isinstance(titles, list) and titles else None
    year = year_of_crossref(cr)
    oq = urllib.parse.quote("https://doi.org/" + doi, safe=":/")
    alex, os = fetch("https://api.openalex.org/works/" + oq)
    if isinstance(alex, dict):
        adoi = str(alex.get("doi") or "").lower()
        if adoi and adoi.removeprefix("https://doi.org/") != doi.lower():
            os, alex = "wrong_doi", {}
    else:
        alex = {}
    access = alex.get("open_access") or {}
    loc = alex.get("best_oa_location") or {}
    if not isinstance(access, dict):
        access = {}
    if not isinstance(loc, dict):
        loc = {}
    url = loc.get("landing_page_url") or loc.get("pdf_url")
    if not isinstance(url, str) or not URL_RE.match(url):
        url = None
    return {
        "crossref_lookup": cs,
        "crossref_title": title,  # public bibliographic metadata
        "crossref_publication_year": year,
        "openalex_lookup": os,
        "openalex_oa_record_claim": bool(access.get("is_oa")),
        "openalex_oa_status": access.get("oa_status") or None,
        "openalex_fulltext_index_claim": bool(alex.get("has_fulltext")),
        "openalex_candidate_oa_location_url": url,
        "candidate_oa_link_NOT_tested_for_article_text": True,
        "actual_visitor_taxa_or_effective_pollination_NOT_inferred_from_abstract": True,
    }


def run(author, preflight, queue_out, summary_out, offline=False):
    for p in (queue_out, summary_out):
        require(not Path(p).exists(), "Never overwrite an already-audited output")
    q, counts = source_queue(author, preflight)
    for j, rec in enumerate(q):
        if rec["original_DOI_literal_status"] == "doi_syntax":
            if offline:
                rec["crossref_lookup"] = "not_queried_offline"
                rec["openalex_lookup"] = "not_queried_offline"
                rec["crossref_title"] = None
                rec["openalex_candidate_oa_location_url"] = None
                rec["openalex_oa_record_claim"] = False
            else:
                rec.update(metadata_for(rec["original_bibliographic_identifier"]))
                time.sleep(.15)
        else:
            rec.update({
                "crossref_lookup": "not_a_DOI",
                "openalex_lookup": "not_a_DOI",
                "openalex_oa_record_claim": False,
                "openalex_candidate_oa_location_url": None,
                "crossref_title": None,
            })
        if (j + 1) % 15 == 0:
            print("ORIGINAL_92_VISITOR_PAPER_METADATA_PROGRESS", j+1, "/", len(q),
                  flush=True)
    summary = {
        "schema": "ALAN_92_real_original_papers_bibliographic_oa_gate_v1",
        "status": "PUBLIC_BIBLIOGRAPHIC_LINKS_ONLY_NO_LOCAL_GUILD",
        "published_original_source_gitblob": ORIGINAL_BLOB,
        "true_reproductive_source_prelight_SHA256": PREFLIGHT_SHA,
        "published_studies_with_author_visitor_observation_flag": 92,
        "genuine_experiment_reproductive_outcome_pairs_of_these_papers": 247,
        "distinct_original_experiment_plants_in_these_papers": 110,
        "original_identifier_class_counts": dict(sorted(counts.items())),
        "valid_DOI_literal_papers_attempted": 87 if not offline else 0,
        "crossref_lookup_counts": dict(sorted(Counter(
            z["crossref_lookup"] for z in q).items())),
        "openalex_lookup_counts": dict(sorted(Counter(
            z["openalex_lookup"] for z in q).items())),
        "crossref_title_returned_count": sum(bool(z["crossref_title"]) for z in q),
        "openalex_OA_claimed_count": sum(bool(
            z["openalex_oa_record_claim"]) for z in q),
        "openalex_candidate_oa_location_url_count": sum(bool(
            z["openalex_candidate_oa_location_url"]) for z in q),
        "papers_with_same_site_year_plant_visitor_guild_confirmed": 0,
        "n_usable_effective_moth_daytime_function_comparisons": 0,
        "metadata_source_does_not_provide_observed_guild_or_effective_pollen": True,
        "public_text_curation_needed_next": True,
        "source_already_literature_selected_not_random_global_sample": True,
    }
    for path in (queue_out, summary_out):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(queue_out).write_text(json.dumps({
        "schema": "ALAN_92_real_study_public_bibliography_curation_queue_v1",
        "source_locked_study_reference_no_plant_site_year": True,
        "records": q,
    }, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    Path(summary_out).write_text(json.dumps(
        summary, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    print("ALAN_REAL_92_SOURCE_PAPER_OA_GATE", json.dumps(
        summary, ensure_ascii=False), flush=True)
    return summary


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("author", "preflight", "queue_out", "summary_out"):
        p.add_argument("--" + k, required=True, type=Path)
    p.add_argument("--offline", action="store_true")
    a = p.parse_args()
    run(a.author, a.preflight, a.queue_out, a.summary_out, a.offline)
