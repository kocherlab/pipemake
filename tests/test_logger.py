import os
import pytest
import tempfile

from pipemake.logger import *


@pytest.mark.parametrize(
    "log_filename",
    ["test.log"],
)
def test_logger_startLogger(log_filename):
    test_dir = tempfile.mkdtemp()

    # Assign the log file path
    log_file_path = os.path.join(test_dir, log_filename)

    # Start the logger
    startLogger(log_filename=log_file_path)

    # Log a message
    logging.info("This is a test message")

    # Check if the log file was created
    assert os.path.exists(log_file_path)


@pytest.mark.parametrize(
    "log_filename",
    ["test.log"],
)
def test_logger_logArgDict(log_filename):
    test_dir = tempfile.mkdtemp()

    # Assign the log file path
    test_dict = {"arg1": "value1", "arg2": None, "arg3": "value3"}

    # Assign the log file path
    log_file_path = os.path.join(test_dir, log_filename)

    # Start the logger
    startLogger(log_filename=log_file_path)

    # Log the argument dictionary
    logArgDict(test_dict)

    # Open the log file
    with open(log_file_path, "r") as log_file:
        # Read the log file
        log_text = log_file.read()

        # Check if the log file contains the expected text
        assert "Argument arg1: value1" in log_text
        assert "Argument arg2: None" not in log_text
        assert "Argument arg3: value3" in log_text


def test_logger_startEarlyLogger_uses_log_format(capsys):
    startEarlyLogger()

    logging.warning("early warning")

    # Same layout as the later console output: "LEVEL: message"
    assert capsys.readouterr().err == "WARNING: early warning\n"


@pytest.mark.parametrize(
    "level, expected",
    [
        (logging.WARNING, "WARNING: warn\n"),
        (logging.INFO, "INFO: info\nWARNING: warn\n"),
        (logging.DEBUG, "DEBUG: debug\nINFO: info\nWARNING: warn\n"),
    ],
)
def test_logger_startEarlyLogger_level(capsys, level, expected):
    startEarlyLogger(level=level)

    logging.debug("debug")
    logging.info("info")
    logging.warning("warn")

    assert capsys.readouterr().err == expected
