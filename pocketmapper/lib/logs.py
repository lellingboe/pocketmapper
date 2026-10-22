"""
Package logging: the record format's `stage` filter and the per-run log file.
"""

import logging
import os
from contextlib import contextmanager

from pocketmapper.lib.constants import LOG_FORMAT, PACKAGE_LOGGER
from pocketmapper.lib.files import make_dir


class StageFilter(logging.Filter):
    """
    Supply a missing `stage` attribute from the name of the function that logged the record.

    `LOG_FORMAT` interpolates `%(stage)s`, which is not a stock LogRecord attribute, so without this
    filter a record logged with no `extra={"stage": ...}` fails to format.
    """

    def filter(self, record):
        """
        Default `record.stage` to `record.funcName` and keep the record.

        Args:
            record (logging.LogRecord): The record about to be emitted; modified in place.

        Returns:
            bool: Always True -- nothing is filtered out.
        """
        if not hasattr(record, "stage"):
            record.stage = record.funcName
        return True


def format_handler(handler):
    """
    Give a handler the package log format and the `StageFilter` that format needs.

    Args:
        handler (logging.Handler): The handler to configure; modified in place.

    Returns:
        logging.Handler: The same handler.
    """
    handler.setFormatter(logging.Formatter(LOG_FORMAT))
    handler.addFilter(StageFilter())
    return handler


@contextmanager
def log_to_file(log_path, verbosity):
    """
    Log the package to a file, at a level set by `verbosity`, for the length of a `with` block.

    Creates the file's directory, sets the `pocketmapper` logger's level and adds a file handler
    appending to `log_path`. Both are undone on exit, however the block ends. If the logger already
    writes to `log_path`, as inside an enclosing call, or `log_path` is None, only the level is set.

    Args:
        log_path (str or None): The log file, or None for none.
        verbosity (int): 4=DEBUG, 3=INFO, 2=WARNING, else ERROR.

    Yields:
        None

    Raises:
        PocketMapperError: If the log's directory cannot be created.
    """
    if log_path is not None:
        make_dir(os.path.dirname(log_path), {"stage": "Configuring Settings"})
    if verbosity == 4:
        log_level = logging.DEBUG
    elif verbosity == 3:
        log_level = logging.INFO
    elif verbosity == 2:
        log_level = logging.WARNING
    else:
        log_level = logging.ERROR

    package_logger = logging.getLogger(PACKAGE_LOGGER)
    previous_log_level = package_logger.level
    package_logger.setLevel(log_level)
    # A second handler on the same file would write every line twice
    handler = None
    if log_path is not None and not any(
        isinstance(existing, logging.FileHandler) and existing.baseFilename == os.path.abspath(log_path)
        for existing in package_logger.handlers
    ):
        handler = format_handler(logging.FileHandler(log_path))
        package_logger.addHandler(handler)
    try:
        yield
    finally:
        if handler is not None:
            package_logger.removeHandler(handler)
            handler.close()
        package_logger.setLevel(previous_log_level)
