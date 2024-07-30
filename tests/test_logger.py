import os
import pytest
import tempfile

from pipemake.logger import *

@pytest.mark.parametrize(
    "log_filename",
    [
        "test.log"
    ],
)
def test_logger_startLogger (log_filename):
    test_dir = tempfile.mkdtemp()

    # Assign the log file path
    log_file_path = os.path.join(test_dir, log_filename)

    # Start the logger
    startLogger(log_filename = log_file_path)

    # Log a message
    logging.info('This is a test message')

    # Check if the log file was created
    assert os.path.exists(log_file_path)

@pytest.mark.parametrize(
    "log_filename",
    [
        "test.log"
    ],
)
def test_logger_logArgDict (log_filename):
    test_dir = tempfile.mkdtemp()

    # Assign the log file path
    test_dict = {'arg1': 'value1', 'arg2': None, 'arg3': 'value3'}

    # Assign the log file path
    log_file_path = os.path.join(test_dir, log_filename)

    # Start the logger
    startLogger(log_filename = log_file_path)

    # Log the argument dictionary
    logArgDict(test_dict)

    # Open the log file
    with open(log_file_path, 'r') as log_file:

        # Read the log file
        log_text = log_file.read()

        # Check if the log file contains the expected text
        assert 'Argument arg1: value1' in log_text
        assert 'Argument arg2: None' not in log_text
        assert 'Argument arg3: value3' in log_text

