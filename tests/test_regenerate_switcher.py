# Copyright (c) 2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-ANL-Open-Source-License
"""Tests for tools/regenerate_switcher.py.

The script lives in ``tools/`` (outside the importable package) so it is
not part of the wheel and not subject to the package coverage threshold.
These tests exercise its three pure functions plus the end-to-end CLI
entry point against a fake gh-pages working tree built in ``tmp_path``.

The script is the source of truth for the published version-switcher
dropdown (issue #273).  Regression here means the documentation site
loses its ability to navigate between published versions, so the
behavior is locked down with explicit tests.
"""

from __future__ import annotations

import importlib.util
import json
import re
import sys
from contextlib import nullcontext as does_not_raise
from pathlib import Path

import pytest

# Load tools/regenerate_switcher.py as a module under a synthetic name.
# It is not on sys.path because tools/ is not a Python package.
_ROOT = Path(__file__).resolve().parents[1]
_SCRIPT = _ROOT / "tools" / "regenerate_switcher.py"
_spec = importlib.util.spec_from_file_location("regenerate_switcher", _SCRIPT)
regen = importlib.util.module_from_spec(_spec)
sys.modules["regenerate_switcher"] = regen
_spec.loader.exec_module(regen)  # type: ignore[union-attr]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_fake_gh_pages(root: Path, names: list[str]) -> None:
    """Create *root*/<name>/_static/ for each entry in *names*."""
    for name in names:
        (root / name / "_static").mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# _is_stable
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name, expected, context",
    [
        pytest.param("v0.10.1", True, does_not_raise(), id="stable-final"),
        pytest.param("v1.0.0", True, does_not_raise(), id="stable-major"),
        pytest.param("v0.4.1rc1", False, does_not_raise(), id="rc-prerelease"),
        pytest.param("v0.4.1a2", False, does_not_raise(), id="alpha-prerelease"),
        pytest.param("v0.4.1b1", False, does_not_raise(), id="beta-prerelease"),
        pytest.param("latest", False, does_not_raise(), id="non-version-name"),
        pytest.param("v0.4", False, does_not_raise(), id="two-component-rejected"),
    ],
)
def test_is_stable(name, expected, context):
    """``_is_stable`` recognizes only well-formed final-release names."""
    with context:
        assert regen._is_stable(name) is expected


# ---------------------------------------------------------------------------
# _version_sort_key (via discover_versions ordering)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "names, expected_order, context",
    [
        pytest.param(
            ["v0.3.0", "v0.10.0", "v0.5.0", "v0.4.1"],
            ["v0.10.0", "v0.5.0", "v0.4.1", "v0.3.0"],
            does_not_raise(),
            id="numeric-not-lexicographic",
        ),
        pytest.param(
            ["v0.4.1", "v0.4.1rc1", "v0.4.1a1"],
            ["v0.4.1", "v0.4.1rc1", "v0.4.1a1"],
            does_not_raise(),
            id="stable-above-prereleases-of-same-base",
        ),
        pytest.param(
            ["v0.10.2rc1", "v0.10.1"],
            ["v0.10.2rc1", "v0.10.1"],
            does_not_raise(),
            id="newer-base-prerelease-above-older-stable",
        ),
    ],
)
def test_discover_versions_sort_order(tmp_path, names, expected_order, context):
    """Versions are returned newest-first with PEP-440-aware ordering."""
    with context:
        _make_fake_gh_pages(tmp_path, names + ["latest"])
        assert regen.discover_versions(tmp_path) == expected_order


@pytest.mark.parametrize(
    "tree, context",
    [
        pytest.param(
            ["latest", "v0.5.0", "not-a-version", "_internal"],
            does_not_raise(),
            id="ignores-non-version-dirs",
        ),
    ],
)
def test_discover_versions_filters_non_version_dirs(tmp_path, tree, context):
    """Directories that don't match ``vX.Y.Z[suffix]`` are ignored."""
    with context:
        _make_fake_gh_pages(tmp_path, tree)
        result = regen.discover_versions(tmp_path)
        assert result == ["v0.5.0"]


def test_discover_versions_missing_dir(tmp_path):
    """Pointing at a non-existent path raises ``FileNotFoundError``."""
    with pytest.raises(FileNotFoundError, match=re.escape("gh-pages directory")):
        regen.discover_versions(tmp_path / "does-not-exist")


# ---------------------------------------------------------------------------
# build_switcher
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "versions, base_url, expected_preferred, context",
    [
        pytest.param(
            ["v0.10.1", "v0.10.0", "v0.9.0"],
            "https://example.com/proj",
            "v0.10.1",
            does_not_raise(),
            id="highest-stable-is-preferred",
        ),
        pytest.param(
            ["v0.10.2rc1", "v0.10.1", "v0.9.0"],
            "https://example.com/proj",
            "v0.10.1",
            does_not_raise(),
            id="prerelease-skipped-for-preferred",
        ),
        pytest.param(
            ["v0.5.0rc1"],
            "https://example.com/proj",
            None,
            does_not_raise(),
            id="no-stable-no-preferred",
        ),
    ],
)
def test_build_switcher_marks_preferred(
    versions, base_url, expected_preferred, context
):
    """The first stable release becomes ``preferred``; pre-releases never do."""
    with context:
        payload = regen.build_switcher(versions, base_url)
        # First entry is always "latest".
        assert payload[0] == {
            "version": "latest",
            "url": f"{base_url}/latest/",
        }
        preferred = [e for e in payload if e.get("preferred")]
        if expected_preferred is None:
            assert preferred == []
        else:
            assert len(preferred) == 1
            assert preferred[0]["version"] == expected_preferred
        # Trailing slash discipline.
        assert all(e["url"].endswith("/") for e in payload)
        # No double slashes from a stray trailing slash on base_url.
        payload2 = regen.build_switcher(versions, base_url + "/")
        assert payload2 == payload


# ---------------------------------------------------------------------------
# write_switcher_files + end-to-end CLI
# ---------------------------------------------------------------------------


def test_write_switcher_files_mirrors_into_every_version(tmp_path):
    """Every vX.Y.Z/_static/ gets the same canonical switcher.json."""
    names = ["latest", "v0.10.1", "v0.10.0", "v0.9.3"]
    _make_fake_gh_pages(tmp_path, names)
    versions = regen.discover_versions(tmp_path)
    payload = regen.build_switcher(versions, "https://example.com/p")
    written = regen.write_switcher_files(tmp_path, payload, versions)

    # latest + 3 versions = 4 files written.
    assert len(written) == len(names)
    canonical = (tmp_path / "latest/_static/switcher.json").read_text()
    for name in names:
        target = tmp_path / name / "_static" / "switcher.json"
        assert target.is_file()
        assert target.read_text() == canonical


def test_write_switcher_files_skips_dirs_without_static(tmp_path, capsys):
    """A version directory missing its _static/ subdir is skipped (defensive)."""
    _make_fake_gh_pages(tmp_path, ["latest", "v0.10.1"])
    # v0.9.0 exists but has no _static/ — discover_versions would normally
    # not return such a dir, but we exercise the defensive branch directly.
    (tmp_path / "v0.9.0").mkdir()
    payload = regen.build_switcher(["v0.10.1", "v0.9.0"], "https://e/p")
    written = regen.write_switcher_files(tmp_path, payload, ["v0.10.1", "v0.9.0"])
    assert len(written) == 2  # latest and v0.10.1, not v0.9.0
    err = capsys.readouterr().err
    assert "skip v0.9.0" in err


def test_main_dry_run(tmp_path, capsys):
    """``--dry-run`` prints the JSON to stdout and writes nothing."""
    names = ["latest", "v0.10.1", "v0.5.0"]
    _make_fake_gh_pages(tmp_path, names)
    rc = regen.main([str(tmp_path), "--dry-run", "--base-url", "https://e/p"])
    assert rc == 0
    out = capsys.readouterr().out
    payload = json.loads(out)
    assert payload[0]["version"] == "latest"
    assert any(e.get("preferred") for e in payload)
    # Nothing was written.
    for name in names:
        assert not (tmp_path / name / "_static" / "switcher.json").exists()


def test_main_writes_files(tmp_path, capsys):
    """The default invocation writes switcher.json into every version dir."""
    names = ["latest", "v0.10.1", "v0.5.0"]
    _make_fake_gh_pages(tmp_path, names)
    rc = regen.main([str(tmp_path), "--base-url", "https://e/p"])
    assert rc == 0
    canonical = (tmp_path / "latest/_static/switcher.json").read_text()
    payload = json.loads(canonical)
    # Highest stable is preferred.
    preferred = [e for e in payload if e.get("preferred")]
    assert len(preferred) == 1
    assert preferred[0]["version"] == "v0.10.1"
    # Mirrored into v0.5.0 too.
    assert (tmp_path / "v0.5.0/_static/switcher.json").read_text() == canonical


def test_main_no_versions_returns_error(tmp_path, capsys):
    """An empty gh-pages tree exits with code 2 and a stderr message."""
    (tmp_path / "latest").mkdir()  # only 'latest' — no version dirs
    rc = regen.main([str(tmp_path)])
    assert rc == 2
    err = capsys.readouterr().err
    assert "no version directories found" in err
