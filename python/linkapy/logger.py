import logging
from pathlib import Path

from rich.logging import RichHandler


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