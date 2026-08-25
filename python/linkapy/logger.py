import logging
import re
import sys
from importlib import metadata
from pathlib import Path

from rich.logging import RichHandler


def _core_dependency_names() -> list[str]:
    '''
    Get dependencies
    '''
    names = []
    for req in metadata.requires("linkapy") or []:
        spec, _, marker = req.partition(";")
        if "extra" in marker:
            continue
        match = re.match(r"[A-Za-z0-9_.\-]+", spec.strip())
        if match:
            names.append(match.group())
    return names


def log_versions(logger: logging.Logger) -> None:
    '''
    debug dump python version and dependency versions
    '''
    logger.debug(f"python: {sys.version.split()[0]}")
    logger.debug(f"linkapy: {metadata.version('linkapy')}")
    for pkg in _core_dependency_names():
        try:
            logger.debug(f"{pkg}: {metadata.version(pkg)}")
        except metadata.PackageNotFoundError:
            logger.warning(f"{pkg}: not installed")


def setup_logger(logfile: Path, verbose: bool = False) -> logging.Logger:
    _logger = logging.getLogger()
    rich_handler = RichHandler(rich_tracebacks=True, show_time=False, show_level=True, show_path=False)
    _fmt = logging.Formatter('%(levelname)s - %(asctime)s - %(message)s', datefmt="%H:%M:%S")
    file_handler = logging.FileHandler(logfile)

    file_handler.setFormatter(_fmt)
    # Set verbosity
    if verbose:
        _logger.setLevel(logging.DEBUG)
        file_handler.setLevel(logging.DEBUG)
        rich_handler.setLevel(logging.DEBUG)
    else:
        _logger.setLevel(logging.INFO)
        file_handler.setLevel(logging.INFO)
        rich_handler.setLevel(logging.INFO)


    _logger.addHandler(rich_handler)
    _logger.addHandler(file_handler)
    return _logger
