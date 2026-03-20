"""Failure logging utilities."""

import json
import traceback
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional


class FailureLogger:
    """Append-only JSONL logger for pipeline failures."""

    def __init__(self, log_dir: Path) -> None:
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.log_path = self.log_dir / "failures.jsonl"

    def log(
        self,
        step: str,
        error: BaseException,
        donor_id: Optional[str] = None,
        celltype: Optional[str] = None,
        window: Optional[int] = None,
        **context: Any,
    ) -> None:
        entry = {
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "step": step,
            "donor_id": donor_id,
            "celltype": celltype,
            "window": window,
            "error_type": type(error).__qualname__,
            "error_msg": str(error),
            "traceback": traceback.format_exception(type(error), error, error.__traceback__),
            "context": context if context else None,
        }
        with open(self.log_path, "a") as f:
            f.write(json.dumps(entry, default=str) + "\n")


@contextmanager
def log_failures(logger: FailureLogger, step: str, **context: Any):
    """Context manager that catches exceptions, logs them, and re-raises."""
    try:
        yield
    except Exception as exc:
        logger.log(step, exc, **context)
        raise


@contextmanager
def skip_and_log(logger: FailureLogger, step: str, **context: Any):
    """Context manager that catches exceptions, logs them, but does NOT re-raise."""
    try:
        yield
    except Exception as exc:
        logger.log(step, exc, **context)
