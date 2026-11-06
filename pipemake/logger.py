import logging
import sys

LOG_FORMAT = "%(asctime)s - %(levelname)s: %(message)s"
LOG_DATEFMT = "%Y-%m-%d %H:%M:%S"
CONSOLE_FORMAT = "%(levelname)s: %(message)s"


def _logUncaught(etype, val, tb):
    # Log the error
    logging.error("%s" % (val), exc_info=(etype, val, tb))


def startEarlyLogger(level=logging.WARNING):
    # Log to the console, in the console format, until the log file exists.
    # Only records at or above level are shown (default: warnings and errors). startLogger() replaces this handler.
    for handler in logging.root.handlers[:]:
        handler.close()
        logging.root.removeHandler(handler)

    early_log = logging.StreamHandler()
    early_log.setLevel(level)
    early_log.setFormatter(logging.Formatter(CONSOLE_FORMAT))

    logging.root.setLevel(min(logging.INFO, level))
    logging.root.addHandler(early_log)

    sys.excepthook = _logUncaught


def startLogger(log_filename=None, filemode="w"):
    # Close any old loggers
    for handler in logging.root.handlers[:]:
        handler.close()
        logging.root.removeHandler(handler)

    # Star the logger
    if log_filename:
        logging.basicConfig(
            filename=log_filename,
            filemode=filemode,
            level="INFO",
            format=LOG_FORMAT,
            datefmt=LOG_DATEFMT,
        )
    else:
        logging.basicConfig(
            stream=sys.stdout,
            level="INFO",
            format=LOG_FORMAT,
            datefmt=LOG_DATEFMT,
        )

    # Start logging to stdout
    stdout_log = logging.StreamHandler()

    # Assign the stdout logging level
    stdout_log.setLevel(logging.WARNING)

    # Define the stdout format
    console_format = logging.Formatter(CONSOLE_FORMAT)

    # Assign the format
    stdout_log.setFormatter(console_format)

    # Add the stdout handler
    logging.getLogger("").addHandler(stdout_log)

    # Update the exception hook to log the exception
    sys.excepthook = _logUncaught


def logArgDict(arg_dict, print_undefined=False, omit=[]):
    # Loop the arguments
    for arg, value in arg_dict.items():
        # Skip arg if in omit list
        if arg in omit:
            continue

        # Report only defined arguments, unless print_undefined is True
        if value is None and not print_undefined:
            continue

        # Log the argument
        logging.info("Argument %s: %s" % (arg, value))
