import argparse
import pytest

from pipemake.parser import PipelineParser, ValueRange, PositiveFloat


def buildArgParser(range_dict, num_type=int):
    """Build an argparse parser from a config-style range definition"""
    arg_args = PipelineParser._processArgRange(
        {"type": num_type, "help": "test arg", "range": range_dict}
    )
    parser = argparse.ArgumentParser()
    parser.add_argument("--value", **arg_args)
    return parser, arg_args


@pytest.mark.parametrize(
    "num_type, min_value, max_value, value, expected",
    [
        (int, 1, 10, "5", 5),
        (int, 1, 10, "1", 1),
        (int, 1, 10, "10", 10),
        (float, 0.0, 1.0, "0.5", 0.5),
        (float, 0.0, 1.0, "1", 1.0),
        (int, 1, None, "1000000", 1000000),
        (int, None, 10, "-100", -100),
        (float, 0.0, None, "0", 0.0),
    ],
)
def test_ValueRange_valid(num_type, min_value, max_value, value, expected):
    result = ValueRange(num_type, min_value, max_value)(value)
    assert result == expected
    assert isinstance(result, num_type)


@pytest.mark.parametrize(
    "num_type, min_value, max_value, value",
    [
        (int, 1, 10, "0"),
        (int, 1, 10, "11"),
        (int, 1, None, "0"),
        (int, None, 10, "11"),
        (float, 0.0, 1.0, "1.5"),
        (float, 0.0, 1.0, "-0.1"),
        (float, 0.0, 1.0, "nan"),
        (float, 0.0, None, "nan"),
        (float, None, 1.0, "nan"),
    ],
)
def test_ValueRange_out_of_range(num_type, min_value, max_value, value):
    with pytest.raises(argparse.ArgumentTypeError):
        ValueRange(num_type, min_value, max_value)(value)


@pytest.mark.parametrize(
    "num_type, value", [(int, "2.5"), (int, "abc"), (float, "abc")]
)
def test_ValueRange_invalid_number(num_type, value):
    with pytest.raises(argparse.ArgumentTypeError, match="not a valid"):
        ValueRange(num_type, 0, 10)(value)


def test_ValueRange_invalid_definition():
    with pytest.raises(TypeError):
        ValueRange(str, 0, 10)
    with pytest.raises(ValueError):
        ValueRange(int)


@pytest.mark.parametrize(
    "range_dict, num_type, valid, invalid, help_text",
    [
        ({"min": 1, "max": 10}, int, "5", "11", "supported range: 1 to 10"),
        ({"min": 1}, int, "1000", "0", "supported range: at least 1"),
        ({"max": 100}, int, "-5", "101", "supported range: at most 100"),
        ({"min": 0.0, "max": 0.5}, float, "0.25", "0.75", "supported range: 0.0 to 0.5"),
        ({"min": 0.0}, float, "3.5", "-0.5", "supported range: at least 0.0"),
    ],
)
def test_processArgRange(range_dict, num_type, valid, invalid, help_text):
    parser, arg_args = buildArgParser(range_dict, num_type)

    assert help_text in arg_args["help"]
    assert "range" not in arg_args
    assert parser.parse_args(["--value", valid]).value == num_type(valid)
    with pytest.raises(SystemExit):
        parser.parse_args(["--value", invalid])


@pytest.mark.parametrize(
    "range_dict",
    [
        {},
        {"min": None},
        {"min": None, "max": None},
        {"min": "one"},
        {"max": "ten"},
        {"min": True},
        {"min": 10, "max": 1},
        {"minimum": 1},
        {"min": 1, "maximum": 10},
    ],
)
def test_processArgRange_invalid_range(range_dict):
    with pytest.raises(Exception):
        buildArgParser(range_dict)


def test_processArgRange_non_dict_range_skipped():
    arg_args = {"type": int, "help": "test arg", "range": "1-10"}
    assert PipelineParser._processArgRange(arg_args) is arg_args


@pytest.mark.parametrize(
    "value, expected",
    [("1", 1.0), ("0.5", 0.5), ("0.001", 0.001), ("2.5", 2.5), ("10", 10.0)],
)
def test_PositiveFloat_valid(value, expected):
    assert PositiveFloat(value) == expected


@pytest.mark.parametrize(
    "value", ["0", "0.0", "-1", "-0.5", "nan", "inf", "-inf", "abc", ""]
)
def test_PositiveFloat_invalid(value):
    with pytest.raises(argparse.ArgumentTypeError):
        PositiveFloat(value)
