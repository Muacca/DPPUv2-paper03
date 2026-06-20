"""
Audit Utilities
===============

Lightweight check-result accumulator for audit-style scripts.

Pairs with ``dppu.utils.tee_logger`` (which handles log-file capture via
stdout tee).  This module is responsible only for *recording* check outcomes
and *writing* them to CSV — it holds no file handles of its own.

Typical usage pattern::

    from dppu.utils.audit import CheckResults
    from dppu.utils.tee_logger import setup_log, teardown_log

    results = CheckResults()

    def main() -> None:
        results.check("C-01", "input file exists", os.path.exists(SOME_PATH))
        results.warn("W-01", "provisional impl", "revisit later")
        # ...
        results.write_csv(CHECK_CSV)

    if __name__ == "__main__":
        setup_log(__file__, log_dir=LOG_DIR)
        try:
            main()
        finally:
            teardown_log()

Default CSV columns: ``check_id``, ``label``, ``status``, ``detail``.
Pass a custom *fields* list to ``__init__`` when the script uses a different
column schema.
"""

from __future__ import annotations

from typing import Sequence

from .io import write_csv

_DEFAULT_FIELDS: list[str] = ["check_id", "label", "status", "detail"]


class CheckResults:
    """Accumulate pass/fail/warn records and write them to CSV.

    Parameters
    ----------
    fields:
        CSV column names in order.  Defaults to
        ``["check_id", "label", "status", "detail"]``.
    """

    def __init__(self, fields: Sequence[str] | None = None) -> None:
        self.fields: list[str] = list(fields) if fields is not None else _DEFAULT_FIELDS[:]
        self.rows: list[dict[str, str]] = []

    # ── core recording ────────────────────────────────────────────────────────

    def record(self, check_id: str, label: str, status: str, detail: str = "") -> None:
        """Append one row and echo it to stdout (captured by tee_logger)."""
        print(f"  [{status}] {check_id}: {label}" + (f" -- {detail}" if detail else ""))
        self.rows.append(
            {
                "check_id": check_id,
                "label": label,
                "status": status,
                "detail": detail,
            }
        )

    def check(self, check_id: str, label: str, passed: bool, detail: str = "") -> bool:
        """Record PASS when *passed* is True, FAIL otherwise.  Returns *passed*."""
        self.record(check_id, label, "PASS" if passed else "FAIL", detail)
        return passed

    def warn(self, check_id: str, label: str, detail: str = "") -> None:
        """Record a WARN entry."""
        self.record(check_id, label, "WARN", detail)

    def check_status(
        self, check_id: str, status: str, detail: str, evidence: str = ""
    ) -> None:
        """Record one row using an explicit *status* string (PASS/FAIL/WARN/…).

        This mirrors the ``Audit.check(cid, status, detail, evidence)`` call
        pattern used in the legacy per-script ``class Audit`` implementations.
        Output goes to stdout (captured by tee_logger when active).
        """
        print(f"  [{status}] {check_id}: {detail}")
        if evidence:
            print(f"    evidence: {evidence}")
        self.rows.append(
            {
                "check_id": check_id,
                "status": status,
                "detail": detail,
                "evidence": evidence,
            }
        )

    def log(self, text: str = "") -> None:
        """Print *text* to stdout (captured by tee_logger when active).

        Drop-in replacement for ``Audit.log(text)`` / ``AuditLog.log(text)``
        in legacy per-script Audit classes.
        """
        print(text)

    # ── output ────────────────────────────────────────────────────────────────

    def write_csv(self, path: str) -> None:
        """Write accumulated rows to *path* using the configured *fields*.

        Parent directories are created automatically (delegated to
        ``dppu.utils.io.write_csv``).
        """
        write_csv(path, self.fields, self.rows)

    # ── convenience ───────────────────────────────────────────────────────────

    def counts(self) -> dict[str, int]:
        """Return a ``{status: count}`` summary dict."""
        result: dict[str, int] = {}
        for row in self.rows:
            s = row.get("status", "")
            result[s] = result.get(s, 0) + 1
        return result

    def has_failures(self) -> bool:
        """Return True if any FAIL record exists."""
        return any(row.get("status") == "FAIL" for row in self.rows)
