"""Tests for ``Search.PRESETS`` and ``Search(preset=...)``.

No database is opened: the behavioural tests hand ``Search`` stub clients, one
of which knows aspirin, so that corroboration (and therefore ``"strict"``) can
be observed without the offline datasets.
"""

import inspect

import pandas as pd
import pytest

from provesid.search import Search

_ASPIRIN_SMILES = "CC(=O)OC1=CC=CC=C1C(=O)O"
_ASPIRIN_IK = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"


class _Empty:
    """An offline source that holds nothing."""

    def __getattr__(self, method):
        return lambda *args, **kwargs: None


class _OfflineAspirin:
    """An offline PubChemID that knows aspirin by every identifier."""

    def __getattr__(self, method):
        row = {"cid": 2244, "smiles": _ASPIRIN_SMILES, "inchikey": _ASPIRIN_IK,
               "cmpdname": "aspirin", "cas_numbers": ["50-78-2"]}
        return lambda *args, **kwargs: [row] if method.startswith("search_") else row


def _one_source_search(**kwargs):
    """A Search on which only PubChemID holds anything: aspirin, uncorroborated."""
    return Search(
        "cas",
        chebi=_Empty(), comptox=_Empty(), pubchem=_OfflineAspirin(), chembl=_Empty(),
        show_progress=False,
        **kwargs,
    )


# ── The table itself ─────────────────────────────────────────────────────────


def test_presets_are_balanced_strict_and_recall():
    assert set(Search.PRESETS) == {"balanced", "strict", "recall"}


def test_every_preset_names_the_same_keys():
    keys = set(Search.PRESETS["balanced"])
    for name, preset in Search.PRESETS.items():
        assert set(preset) == keys, name


def test_every_preset_key_is_a_constructor_argument_defaulting_to_none():
    params = inspect.signature(Search.__init__).parameters
    for key in Search.PRESETS["balanced"]:
        assert key in params, key
        assert params[key].default is None, key


def test_balanced_is_the_default():
    s = Search("cas")
    assert s.preset == "balanced"
    assert s.settings == Search.PRESETS["balanced"]


def test_strict_differs_from_balanced_only_in_source_support():
    diff = {
        k: v for k, v in Search.PRESETS["strict"].items()
        if Search.PRESETS["balanced"][k] != v
    }
    assert diff == {"min_source_support": 2}


def test_recall_widens_every_way():
    s = Search("name", preset="recall")
    assert s.fuzzy and s.inchikey_skeleton and s.use_zeropm
    assert s.similarity_threshold == 0.7
    assert s.n_hits == "all"
    assert "zeropm" in s._SOURCE_KEYS


@pytest.mark.parametrize("name", sorted(Search.PRESETS))
def test_settings_reproduce_the_preset(name):
    assert Search("cas", preset=name).settings == Search.PRESETS[name]


# ── Explicit arguments override ──────────────────────────────────────────────


def test_explicit_argument_overrides_the_preset():
    s = Search("name", preset="strict", n_hits=3)
    assert s.n_hits == 3
    assert s.min_source_support == 2


def test_explicit_argument_equal_to_the_balanced_default_still_overrides():
    # The case a plain keyword default could not tell apart from "not passed".
    s = Search("name", preset="recall", fuzzy=False, use_zeropm=False)
    assert s.fuzzy is False
    assert s.use_zeropm is False
    assert "zeropm" not in s._SOURCE_KEYS
    assert s.inchikey_skeleton is True


def test_explicit_zero_overrides_strict_support():
    assert Search("cas", preset="strict", min_source_support=0).min_source_support == 0


def test_unknown_preset_is_refused():
    with pytest.raises(ValueError, match="preset must be one of"):
        Search("cas", preset="lenient")


def test_invalid_scorer_is_still_refused_under_a_preset():
    with pytest.raises(ValueError, match="fuzzy_scorer"):
        Search("name", preset="recall", fuzzy_scorer="nonsense")


def test_presets_are_not_aliased_by_an_instance():
    s = Search("cas", preset="strict")
    s.settings["min_source_support"] = 99
    s.min_source_support = 5
    assert Search.PRESETS["strict"]["min_source_support"] == 2
    assert s.settings["min_source_support"] == 5


# ── What the presets do, and what the frame records ─────────────────────────


def test_strict_drops_an_uncorroborated_hit_balanced_keeps():
    balanced = _one_source_search().search("50-78-2")
    strict = _one_source_search(preset="strict").search("50-78-2")
    assert balanced.loc[0, "InChIKey"] == _ASPIRIN_IK
    assert pd.isna(strict.loc[0, "InChIKey"])


def test_attrs_record_the_preset_and_the_call_settings():
    s = _one_source_search(preset="strict")
    df = s.search("50-78-2", min_source_support=1)
    assert df.attrs["preset"] == "strict"
    assert df.attrs["settings"]["min_source_support"] == 1
    assert df.attrs["settings"]["n_hits"] == 1
    assert df.loc[0, "InChIKey"] == _ASPIRIN_IK
    # The per-call override does not stick to the instance.
    assert s.settings["min_source_support"] == 2


def test_enrich_carries_the_preset_and_settings():
    s = _one_source_search(preset="strict")
    out = s.enrich(pd.DataFrame({"CAS": ["50-78-2"]}), "CAS", n_hits=2)
    assert out.attrs["preset"] == "strict"
    assert out.attrs["settings"]["n_hits"] == 2
    assert out.attrs["settings"]["min_source_support"] == 2
