#!/usr/bin/env python3
"""regenerate_switcher.py — rebuild switcher.json from a gh-pages tree.

The pydata-sphinx-theme version-switcher dropdown is driven by the JSON
file at ``latest/_static/switcher.json``.  Historically, that file was
maintained by an *incremental* tag-push step that cloned gh-pages,
inserted the freshly built version, and pushed back.  When that step
failed (or was skipped) for any release, the live file silently fell
behind — which is exactly what happened between v0.5.0 and v0.10.1.

This script replaces the incremental approach with a regeneration
strategy:  the **directory listing of the gh-pages branch is the
source of truth.**  Whatever ``vX.Y.Z`` directories exist on gh-pages
end up in the dropdown, sorted newest-first, with the highest stable
release marked as ``preferred``.  Pre-release versions (``aN``,
``bN``, ``rcN``) appear in the dropdown but never become preferred.

The generated JSON is written to ``<gh-pages>/latest/_static/switcher.json``
and *mirrored* into every ``<gh-pages>/vX.Y.Z/_static/switcher.json`` so
that visitors browsing an older version see the same up-to-date dropdown.

Usage::

    python3 tools/regenerate_switcher.py <gh-pages-checkout> \\
        --base-url https://BCDA-APS.github.io/ad_hoc_diffractometer

The script is idempotent and safe to re-run.  It exits non-zero on any
filesystem error or if no version directories are found.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

# A version directory name like "v0.10.1", "v0.4.1rc1", "v1.2.3a2".
# Capture groups: major, minor, patch, optional pre-release tag (a|b|rc + N).
_VERSION_RE = re.compile(r"^v(\d+)\.(\d+)\.(\d+)((?:a|b|rc)\d+)?$")


def _version_sort_key(name: str) -> tuple:
    """Return a tuple suitable for sorting version directory names.

    Stable releases sort *after* their pre-releases (so v0.4.1 ranks
    above v0.4.1rc1 when sorting newest-first via ``reverse=True``).
    The tuple shape ``(major, minor, patch, stable_flag, pre_kind, pre_n)``
    preserves the natural numeric order; ``stable_flag`` is 1 for a
    final release and 0 for a pre-release, so a stable release of the
    same base version sorts higher.
    """
    m = _VERSION_RE.match(name)
    if not m:
        # Should never happen — caller filters by the same regex.
        return (-1, -1, -1, -1, "", -1)
    major, minor, patch = int(m.group(1)), int(m.group(2)), int(m.group(3))
    pre = m.group(4)
    if pre is None:
        return (major, minor, patch, 1, "", 0)
    kind = pre.rstrip("0123456789")  # "a", "b", or "rc"
    n = int(pre[len(kind) :])
    # a < b < rc lexicographically already; that matches PEP 440 ordering.
    return (major, minor, patch, 0, kind, n)


def _is_stable(name: str) -> bool:
    """Return True if *name* is a final release (no a/b/rc suffix)."""
    m = _VERSION_RE.match(name)
    return m is not None and m.group(4) is None


def discover_versions(gh_pages: Path) -> list[str]:
    """Return the sorted (newest-first) list of vX.Y.Z dirs on gh-pages."""
    if not gh_pages.is_dir():
        raise FileNotFoundError(f"gh-pages directory not found: {gh_pages}")
    names = [
        p.name for p in gh_pages.iterdir() if p.is_dir() and _VERSION_RE.match(p.name)
    ]
    names.sort(key=_version_sort_key, reverse=True)
    return names


def build_switcher(versions: list[str], base_url: str) -> list[dict]:
    """Build the switcher.json payload as a list of dicts.

    The first entry is always ``latest``; remaining entries are the
    version directories newest-first.  The highest **stable** version
    is marked ``preferred`` so the pydata-sphinx-theme banner that
    points users at the stable release links to the right place.
    """
    base_url = base_url.rstrip("/")
    entries: list[dict] = [
        {"version": "latest", "url": f"{base_url}/latest/"},
    ]

    preferred_marked = False
    for v in versions:
        entry: dict = {"version": v, "url": f"{base_url}/{v}/"}
        if not preferred_marked and _is_stable(v):
            entry["preferred"] = True
            preferred_marked = True
        entries.append(entry)

    return entries


def write_switcher_files(
    gh_pages: Path, payload: list[dict], versions: list[str]
) -> list[Path]:
    """Write switcher.json into latest/ and every version directory.

    Returns the list of files written.  Each version dir gets the
    *same* canonical payload so the dropdown looks identical no matter
    which version page the visitor is viewing.
    """
    text = json.dumps(payload, indent=4) + "\n"
    written: list[Path] = []

    targets = ["latest"] + versions
    for name in targets:
        static_dir = gh_pages / name / "_static"
        if not static_dir.is_dir():
            # No _static dir means this isn't a Sphinx build output;
            # skip silently rather than fail (defensive — should not
            # happen for any directory matched by _VERSION_RE).
            print(f"  skip {name}: no _static directory", file=sys.stderr)
            continue
        target = static_dir / "switcher.json"
        target.write_text(text, encoding="utf-8")
        written.append(target)
    return written


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "gh_pages",
        type=Path,
        help="Path to a checked-out gh-pages working tree.",
    )
    p.add_argument(
        "--base-url",
        default="https://BCDA-APS.github.io/ad_hoc_diffractometer",
        help="Base URL where docs are served (default: %(default)s).",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the generated switcher.json to stdout without writing.",
    )
    args = p.parse_args(argv)

    versions = discover_versions(args.gh_pages)
    if not versions:
        print(
            f"error: no version directories found under {args.gh_pages}",
            file=sys.stderr,
        )
        return 2

    payload = build_switcher(versions, args.base_url)

    if args.dry_run:
        print(json.dumps(payload, indent=4))
        return 0

    written = write_switcher_files(args.gh_pages, payload, versions)
    print(
        f"Regenerated switcher.json from {len(versions)} version "
        f"directories; wrote {len(written)} file(s)."
    )
    for path in written:
        print(f"  {path.relative_to(args.gh_pages)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
