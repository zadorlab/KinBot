"""Module for the configuration of how and what is recorded in the log file."""
import sys
import os
import logging
import warnings
from importlib.metadata import version
from pathlib import Path
import tomllib
from dateutil import parser

from kinbot import license_message
from kinbot.rdkit_config import log_rdkit


def _kinbot_version():
    """Use the running checkout's version, or installed package metadata."""
    project = Path(__file__).resolve().parents[1] / 'pyproject.toml'
    if project.is_file():
        with project.open('rb') as stream:
            return tomllib.load(stream)['project']['version']
    return version('kinbot')


def log_exception(exc_type, exc_value, exc_tb):
    """Sets up the recording of exceptions on the log file

    @param exc_type: Type of exception
    @param exc_value: Value of the exception
    @param exc_tb:
    @return: None
    """
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_tb)
        return
    logger = logging.getLogger('KinBot')
    logger.error("", exc_info=(exc_type, exc_value, exc_tb))


def log_warning(message, *args, **kwargs):
    """Sets up the recording of warnings on the log file

    @param message: Warning message.
    @param args: Additional arguments.
    @param kwargs: Additional keyword arguments.
    @return: None
    """
    logger = logging.getLogger('KinBot')
    logger.warning(" ".join(f"{message}".split()))


def config_log(label, mode='kinbot', level='info'):
    # TODO Format log to break long lines (after column 80).
    """Configures the logger to record all calculation events on a log file.

    @param label: Label of the logger to be used.
    @param mode: The mode of the logger.
    @param level: the level of verbosity for the log file.
    @return: The logger object.
    """
    logger = logging.getLogger(label)
    if level == 'info':
        logger.setLevel(logging.INFO)
    elif level == 'debug':
        logger.setLevel(logging.DEBUG)
    fname = f'{mode}.log'

    # Verbose setup configures this logger again. Do not keep writing the
    # current run's header and messages into its archived log.
    for handler in list(logger.handlers):
        if getattr(handler, '_kinbot_run_log', False):
            logger.removeHandler(handler)
            handler.close()

    # Backup previous log
    if os.path.isfile(fname):
        with open(fname) as log_fh:
            first_line = log_fh.readline()
        try:
            dt = parser.parse(first_line.partition('-INFO:')[0])
        except ValueError:
            dt = 'old'
        date = str(dt).replace(' ', '_')
        os.rename(fname, f'{mode}_{date}.log')
    log_handler = logging.FileHandler(fname, mode='w')
    log_handler._kinbot_run_log = True
    if level in ['debug', 'verbose']:
        log_handler.setLevel(logging.DEBUG)
    else:
        log_handler.setLevel(logging.INFO)

    log_format = logging.Formatter(fmt='%(asctime)s-%(levelname)s: %(message)s',
                                   datefmt='%d-%b-%y %H:%M:%S')
    log_handler.setFormatter(log_format)

    logger.addHandler(log_handler)

    sys.excepthook = log_exception
    warnings.showwarning = log_warning
    logger.info('KinBot %s', _kinbot_version())
    log_rdkit(logger)
    logger.info(license_message.message)

    return logger
