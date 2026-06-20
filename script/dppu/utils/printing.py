"""
DPPU Printing Utilities
=======================

Console output helpers for consistent formatting across scripts.

These functions are intentionally simple so that standalone scripts
can reproduce the same look without importing this module.

Author: Muacca
"""

import json
from typing import Any, Sequence


def hline(char: str = '=', width: int = 70) -> None:
    """Print a horizontal rule."""
    print(char * width)


def print_header(title: str, width: int = 70) -> None:
    """Print a section header with surrounding rules."""
    print()
    hline(width=width)
    print(f"  {title}")
    hline(width=width)


def print_sub(title: str) -> None:
    """Print a sub-section title."""
    print()
    print(f"  \u2500\u2500 {title} \u2500\u2500")
    print()


# ── Markdown ──────────────────────────────────────────────────────────────────

def md_table(
    rows: "Sequence[dict]",
    fields: "Sequence[str]",
) -> str:
    """Render a list of row dicts as a Markdown table string.

    Args:
        rows:   Sequence of row dicts.  Missing keys are rendered as empty
                cells.  Newlines inside cell values are replaced with
                ``<br>``.
        fields: Column names in order (also used as header labels).

    Returns:
        Multi-line Markdown table string, or ``"_No rows._"`` when *rows*
        is empty.

    Example::

        print(md_table([{"a": 1, "b": 2}], ["a", "b"]))
        # | a | b |
        # | --- | --- |
        # | 1 | 2 |
    """
    if not rows:
        return "_No rows._"
    head = "| " + " | ".join(fields) + " |"
    sep = "| " + " | ".join("---" for _ in fields) + " |"
    body = [
        "| "
        + " | ".join(
            str(row.get(f, "")).replace("\n", "<br>") for f in fields
        )
        + " |"
        for row in rows
    ]
    return "\n".join([head, sep] + body)


# ── Symbolic matrix display ────────────────────────────────────────────────────

def mat_str(mat: Any) -> str:
    """Return a compact string representation of a SymPy matrix.

    Entries are separated by commas within a row and semicolons between
    rows.  Each entry is simplified with ``canonical_simplify`` before
    formatting.  The result is enclosed in square brackets::

        [a,b; c,d]

    Args:
        mat: SymPy Matrix.

    Returns:
        Compact string, e.g. ``"[1,0,0; 0,1,0; 0,0,1]"``.
    """
    from dppu.utils.symbolic import canonical_simplify
    return "[" + "; ".join(
        ",".join(str(canonical_simplify(mat[i, j])) for j in range(mat.cols))
        for i in range(mat.rows)
    ) + "]"


# ── Machine-readable JSON row output ──────────────────────────────────────────

def print_json_rows(
    label: str,
    rows: Sequence[dict],
    fields: Sequence[str],
) -> None:
    """Print rows as machine-readable JSON lines framed by begin/end markers.

    Output protocol::

        <label>_begin
          <label>_row {"field1": "val", ...}
          ...
        <label>_end

    The JSON dict preserves *fields* insertion order.  Values missing from
    a row are rendered as the empty string ``""``.  ``json.dumps`` is called
    with ``sort_keys=False, ensure_ascii=True``.

    Args:
        label:  Prefix token (no spaces).
        rows:   Sequence of row dicts.
        fields: Column names in the desired order.
    """
    print(f"{label}_begin")
    for row in rows:
        ordered = {field: row.get(field, "") for field in fields}
        print(f"  {label}_row " + json.dumps(ordered, sort_keys=False, ensure_ascii=True))
    print(f"{label}_end")