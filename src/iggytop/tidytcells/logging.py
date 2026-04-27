"""Logging configuration for tidytcells warnings (deduplicated file handler)."""

import logging
import os
from datetime import datetime
from logging import getLogger


class DeduplicatingTTFileHandler(logging.FileHandler):
    """File handler that writes each tidytcells warning message only once."""

    def __init__(self, filename: str):
        super().__init__(filename, mode="a", encoding="utf-8")
        self._seen_lines = set()
        if os.path.exists(filename):
            with open(filename, encoding="utf-8") as existing_file:
                for line in existing_file:
                    normalized = line.strip()
                    if normalized:
                        self._seen_lines.add(normalized)

    def emit(self, record: logging.LogRecord) -> None:
        rendered = self.format(record).strip()
        if not rendered:
            return
        self.acquire()
        try:
            if rendered in self._seen_lines:
                return
            self._seen_lines.add(rendered)
            super().emit(record)
        finally:
            self.release()


# Configure tidytcells warning logging on import
log_dir = os.getenv("IGGYTOP_LOG_DIR", "biocypher-log")
os.makedirs(log_dir, exist_ok=True)
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
log_filename = f"tt_standardization_warnings_{timestamp}.log"
log_path = os.path.join(log_dir, log_filename)

tt_logger = getLogger("tidytcells")
tt_logger.setLevel(logging.WARNING)
tt_logger.propagate = False

if not any(isinstance(handler, DeduplicatingTTFileHandler) for handler in tt_logger.handlers):
    handler = DeduplicatingTTFileHandler(log_path)
    handler.setLevel(logging.WARNING)
    handler.setFormatter(logging.Formatter("%(name)s|%(levelname)s|%(message)s"))
    tt_logger.addHandler(handler)
print(f"Find logs related to tidytcells standardization in the {log_path} file.")
