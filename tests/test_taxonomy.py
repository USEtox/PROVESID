"""Tests for the chebifier taxonomy backend (provesid.taxonomy).

Most tests run without the optional ``chebifier`` extra installed: they exercise
the guard/feature-detection, the storage redirect, prediction normalisation, the
tidy-schema row building, and ``to_labels``. A single end-to-end test is skipped
unless ``chebifier`` is actually installed.
"""

import os

import pandas as pd
import pytest

from provesid.taxonomy import (
    CHEBIFIER_PINNED_VERSION,
    TAXONOMY_COLUMNS,
    ChebifierClassifier,
    ChebifierMissingError,
    _configure_chebifier_storage,
    chebifier_available,
    classify_chebifier,
    ensure_v244_indices,
)


class TestFeatureDetection:
    """The optional dependency must be detectable and guarded cleanly."""

    def test_chebifier_available_is_bool(self):
        assert isinstance(chebifier_available(), bool)

    def test_missing_extra_raises_actionable_error(self):
        """Without chebifier, classify raises ChebifierMissingError (not ImportError)."""
        if chebifier_available():
            pytest.skip("chebifier is installed; missing-extra path not exercised")
        # use_cache=False forces the ensemble to load (bypassing the shared,
        # cross-environment on-disk cache), so the guard is actually exercised.
        clf = ChebifierClassifier(use_cache=False)
        with pytest.raises(ChebifierMissingError) as exc:
            clf.classify(["c1ccccc1"])
        assert "install" in str(exc.value).lower()


class TestStorageRedirect:
    """Model caches must be redirected under the PROVESID dataset dir."""

    def test_configure_sets_env_when_unset(self, tmp_path, monkeypatch):
        for var in ("HF_HOME", "HF_HUB_CACHE", "TORCH_HOME"):
            monkeypatch.delenv(var, raising=False)
        base = _configure_chebifier_storage(str(tmp_path))
        assert base == str(tmp_path)
        assert os.environ["HF_HOME"] == os.path.join(str(tmp_path), "huggingface")
        assert os.environ["TORCH_HOME"] == os.path.join(str(tmp_path), "torch")

    def test_configure_respects_existing_env(self, tmp_path, monkeypatch):
        monkeypatch.setenv("HF_HOME", "/already/set")
        _configure_chebifier_storage(str(tmp_path))
        assert os.environ["HF_HOME"] == "/already/set"

    def test_configure_honors_provesid_data_dir(self, tmp_path, monkeypatch):
        for var in ("HF_HOME", "HF_HUB_CACHE", "TORCH_HOME"):
            monkeypatch.delenv(var, raising=False)
        monkeypatch.setenv("PROVESID_DATA_DIR", str(tmp_path))
        base = _configure_chebifier_storage(None)
        assert str(tmp_path) in base


class TestIndexPatch:
    """The v244 index compatibility helper must be safe to call anywhere."""

    def test_ensure_v244_indices_noop_without_chebai_graph(self):
        if pytest.importorskip is None:  # pragma: no cover
            pass
        try:
            import chebai_graph  # noqa: F401
            installed = True
        except ImportError:
            installed = False
        result = ensure_v244_indices()
        if not installed:
            assert result == {}
        else:
            # Idempotent: after one call, a second call reports everything "ok".
            ensure_v244_indices()
            second = ensure_v244_indices()
            assert all(v in {"ok", "missing"} for v in second.values())


class TestPredictionNormalisation:
    """Different chebifier return shapes collapse to {id: confidence}."""

    def test_none_prediction(self):
        assert ChebifierClassifier._normalize_prediction(None) == {}

    def test_list_of_ids(self):
        norm = ChebifierClassifier._normalize_prediction(["22712", 33655])
        assert norm == {"22712": None, "33655": None}

    def test_dict_with_scores(self):
        norm = ChebifierClassifier._normalize_prediction({"22712": 0.9, "33655": None})
        assert norm == {"22712": 0.9, "33655": None}


class TestRowAndLabels:
    """Row assembly and to_labels operate on plain data (no model needed)."""

    def _classifier(self):
        return ChebifierClassifier(use_cache=False, resolve_names=False)

    def test_build_row_schema_and_join(self):
        row = self._classifier()._build_row(
            "c1ccccc1", "UHOVQNZJYSORNB-UHFFFAOYSA-N",
            {"22712": None, "33655": None},
        )
        assert set(row) == set(TAXONOMY_COLUMNS)
        assert row["source"] == "chebifier"
        assert row["chebi_ids"] == "22712|33655"
        assert row["chebi_names"] is None
        assert row["confidence"] is None
        # ClassyFire levels are empty for this backend
        assert row["kingdom"] is None and row["class"] is None

    def test_build_row_with_confidence(self):
        row = self._classifier()._build_row(
            "c1ccccc1", "KEY", {"22712": 0.9, "33655": 0.5}
        )
        assert row["confidence"] == "0.9|0.5"

    def test_build_row_empty_prediction(self):
        row = self._classifier()._build_row("c1ccccc1", "KEY", {})
        assert row["chebi_ids"] is None

    def test_to_labels(self):
        df = pd.DataFrame(
            [{"inchikey": "KEY1", "chebi_ids": "a|b"},
             {"inchikey": "KEY2", "chebi_ids": "c"}]
        )
        labels = ChebifierClassifier.to_labels(df, level="chebi_ids")
        assert labels == {"KEY1": "a|b", "KEY2": "c"}

    def test_to_labels_unknown_level_raises(self):
        df = pd.DataFrame([{"inchikey": "KEY", "chebi_ids": "a"}])
        with pytest.raises(KeyError):
            ChebifierClassifier.to_labels(df, level="nope")


class TestPinnedVersion:
    def test_pinned_version_constant(self):
        assert CHEBIFIER_PINNED_VERSION == "1.2.1"


@pytest.mark.chebifier
@pytest.mark.skipif(
    not chebifier_available(), reason="chebifier extra not installed"
)
class TestLiveClassification:
    """End-to-end classification; only runs when chebifier is installed."""

    def test_classify_benzene_returns_tidy_table(self, tmp_path):
        clf = ChebifierClassifier(data_dir=str(tmp_path), use_cache=True)
        df = clf.classify(["c1ccccc1"])
        assert list(df.columns) == TAXONOMY_COLUMNS
        assert len(df) == 1
        assert df.iloc[0]["source"] == "chebifier"
        assert df.iloc[0]["inchikey"] == "UHOVQNZJYSORNB-UHFFFAOYSA-N"
        assert df.iloc[0]["chebi_ids"]  # non-empty

    def test_convenience_function(self, tmp_path):
        df = classify_chebifier(["c1ccccc1"], data_dir=str(tmp_path))
        assert len(df) == 1
