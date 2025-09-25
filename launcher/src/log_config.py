"""
Logging Module to support proper Python logging.
"""
import sys
import logging
import logging.config


LOGGING = {
    "version": 1,
    "disable_existing_loggers": False,
    "formatters": {
        "verbose": {
            "format": "%(asctime)s [%(filename)s:%(module)s:%(funcName)s:%(lineno)d] %(levelname)s: %(message)s"
        },
        "simple": {"format": "%(levelname)s %(message)s"},
    },
    "handlers": {
        "console": {
            "level": "DEBUG",
            "class": "logging.StreamHandler",
            "formatter": "simple",
        },
        "file": {
            "level": "DEBUG",
            "class": "logging.handlers.TimedRotatingFileHandler",
            "filename": "",
            "when": "midnight",
            "backupCount": 100,
            "formatter": "verbose",
        },
    },
    "loggers": {
        "": {
            "handlers": ["file"],
            "propagate": True,
            "level": "INFO",
        }
    },
}


def enable_logging(log_file="launcher.log", log_level="INFO"):
    """
    Logging configuration.
    """
    LOGGING["handlers"]["file"]["filename"] = log_file

    logging.config.dictConfig(LOGGING)
    _logger = logging.getLogger()
    _logger.setLevel(level=log_level)

    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(level=log_level)

    _logger.addHandler(handler)
    return _logger
