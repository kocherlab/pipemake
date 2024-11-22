import os
import re
import sys
import random
import string
import logging
import argparse
import datetime

from pydoc import locate


def jobRandomString(num_chars=4):
    global random_string

    # Create a random string, if not done already
    if not random_string:
        random_string = "".join(
            random.choices(string.digits + string.ascii_uppercase, k=num_chars)
        )

    return random_string


def jobTimeStamp():
    global time_stamp

    # Create a time stamp, if not done already
    if not time_stamp:
        time_stamp = datetime.datetime.now().strftime("%Y-%m-%d")

    return time_stamp


def confirmDir():
    """Custom action to confirm directory exists"""

    class customAction(argparse.Action):
        def __call__(self, parser, args, value, option_string=None):
            if not os.path.isdir(value):
                raise IOError(f"Unable to find directory: {value}")
            setattr(args, self.dest, value)

    return customAction


def confirmFile():
    """Custom action to confirm file exists"""

    class customAction(argparse.Action):
        def __call__(self, parser, args, value, option_string=None):
            if not os.path.isfile(value):
                raise IOError(f"Unable to find file: {value}")
            setattr(args, self.dest, value)

    return customAction


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write("error: %s\n\n" % message)
        self.print_help()
        sys.exit(2)


class SubcommandHelpFormatter(argparse.RawDescriptionHelpFormatter):
    def _format_action(self, action):
        parts = super(argparse.RawDescriptionHelpFormatter, self)._format_action(action)
        if action.nargs == argparse.PARSER:
            parts = "\n".join(parts.split("\n")[1:])
        return parts


class PipelineParser:
    def __init__(self, config_pipelines):
        # Create the pipeline selection parser
        self._pipeline_parser = MyParser(formatter_class=SubcommandHelpFormatter)
        self._pipeline_parser._positionals.title = (
            "Pipeline selection argument options (positional)"
        )

        # Create the pipeline subparsers
        self._pipeline_subparsers = self._pipeline_parser.add_subparsers(
            dest="pipeline", required=True
        )

        # Ceate a dictionary to store each subparser
        self._pipeline_subparsers_dict = {}

        for pipeline_name, pipeline_config in config_pipelines.items():
            # Assign the arguments for each pipeline
            self._assignPipeline(pipeline_name, pipeline_config)

    def _assignPipeline(self, pipeline_name, pipeline_config):
        # Assign the pipeline parser
        pipeline_parser_dict = pipeline_config.parser_dict

        # Create the subparser
        self._pipeline_subparsers_dict[pipeline_name] = (
            self._pipeline_subparsers.add_parser(
                pipeline_name,
                help=pipeline_parser_dict["help"] + f" (v{pipeline_config.version})",
                add_help=False,
            )
        )

        # Create the argument groups for the pipeline
        pipeline_arg_groups = {}

        # Add the required argument group
        pipeline_arg_groups["required"] = self._pipeline_subparsers_dict[
            pipeline_name
        ].add_argument_group(f"{pipeline_name} required arguments")

        # Loop the pipeline groups
        for pipeline_arg_group in pipeline_parser_dict["arg-groups"]:
            # Create the argument group, if not the basic subparsers
            if pipeline_arg_group == "basic":
                continue

            # Check if the group already exists, then raise an exception if it does
            if pipeline_arg_group in pipeline_arg_groups:
                raise Exception(
                    f"Cannot create group {pipeline_arg_group}. Group already exists."
                )

            # Create the argument group
            pipeline_arg_groups[pipeline_arg_group] = self._pipeline_subparsers_dict[
                pipeline_name
            ].add_argument_group(f"{pipeline_name} {pipeline_arg_group} arguments")

        # Add the optional argument group
        pipeline_arg_groups["optional"] = self._pipeline_subparsers_dict[
            pipeline_name
        ].add_argument_group(f"{pipeline_name} optional arguments")

        # Loop the pipeline argument groups
        for pipeline_arg_group, group_args in pipeline_parser_dict[
            "arg-groups"
        ].items():
            # Check if mutually exclusive arguments are present
            for mutually_exclusive_name, mutually_exclusive_args in self._yeildGroup(
                "mutually-exclusive-args", group_args
            ):
                # Confirm the mutually exclusive group has not been created
                if mutually_exclusive_name in pipeline_arg_groups:
                    raise Exception(
                        f"Cannot create mutually exclusive group {mutually_exclusive_name}. Group already exists."
                    )

                # Check if the group is required, if not given
                if "required" not in mutually_exclusive_args:
                    mutually_exclusive_args["required"] = False

                # Reassign the pipeline_arg_group based on required status
                if (
                    pipeline_arg_group == "basic"
                    and "required" in mutually_exclusive_args
                ):
                    mutually_exclusive_group = "required"
                elif pipeline_arg_group == "basic":
                    mutually_exclusive_group = "optional"
                else:
                    mutually_exclusive_group = pipeline_arg_group

                # Create the mutually exclusive group
                pipeline_arg_groups[mutually_exclusive_name] = pipeline_arg_groups[
                    mutually_exclusive_group
                ].add_mutually_exclusive_group(
                    required=mutually_exclusive_args["required"]
                )

            # Create dict to store the arg wildcard
            arg_wildcards = {}

            # Check if wildcards are present
            for wildcard_name, wildcards_args in self._yeildGroup(
                "wildcards-args", group_args
            ):
                # Process the default arg of the wildcard
                wildcards_args = self._processArgDefaults(wildcards_args)

                # Assign the argument wildcards, if applicable
                arg_wildcards[wildcard_name] = re.findall(
                    r"\{(.*?)\}", wildcards_args["default"]
                )

                # Add the argument to the args argument group
                group_args["args"][wildcard_name] = wildcards_args

            # Check if wildcards are present
            for pipeline_arg_name, arg_args in self._yeildGroup("args", group_args):
                # Process and arg type, if applicable
                if "type" in arg_args:
                    arg_args = self._processArgType(arg_args)

                # Process the arg action, if applicable
                if "action" in arg_args:
                    arg_args = self._processArgAction(arg_args)

                # Process the arg default, if applicable
                if "default" in arg_args:
                    arg_args = self._processArgDefaults(arg_args)

                # Assign the argument wildcards, if applicable
                if "wildcards" in arg_args:
                    arg_args = self._processArgWildcard(
                        pipeline_arg_name, arg_args, arg_wildcards
                    )

                # Check if a default value exists
                elif "default" in arg_args and not isinstance(
                    arg_args["default"], dict
                ):
                    # Add the default string to the help message
                    arg_args["help"] += f' (default: { arg_args["default"]})'

                # Reassign the pipeline_arg_group based on required status
                if "mutually-exclusive" in arg_args:
                    arg_group = arg_args["mutually-exclusive"]
                    del arg_args["mutually-exclusive"]
                elif pipeline_arg_group == "basic" and "required" in arg_args:
                    arg_group = "required"
                elif pipeline_arg_group == "basic":
                    arg_group = "optional"
                else:
                    arg_group = pipeline_arg_group

                # Assign the argument
                pipeline_arg_groups[arg_group].add_argument(
                    f"--{pipeline_arg_name}", **arg_args
                )

        # Add the common optional arguments, but at the end of the list
        pipeline_arg_groups["optional"].add_argument(
            "--workflow-prefix",
            help="Assign the workflow prefix",
            type=str,
            default=f"Workflow_{jobTimeStamp()}_{jobRandomString()}",
        )
        pipeline_arg_groups["optional"].add_argument(
            "--work-dir",
            help="Assign the working directory for snakemake. If not provided, the current directory will be used.",
            type=str,
        )
        pipeline_arg_groups["optional"].add_argument(
            "--scale-threads",
            help="Scale the threads for each task",
            type=float,
            default=1.0,
        )
        pipeline_arg_groups["optional"].add_argument(
            "--scale-mem",
            help="Scale the memory (RAM) for each task",
            type=float,
            default=1.0,
        )
        pipeline_arg_groups["optional"].add_argument(
            "--resource-yml",
            help="Create a seperate resource yaml",
            action="store_true",
        )
        pipeline_arg_groups["optional"].add_argument(
            "--singularity-dir",
            help="Assign different directory of singularity images",
            type=str,
        )
        pipeline_arg_groups["optional"].add_argument(
            "--no-overwrite",
            help="Do not overwrite existing files",
            action="store_false",
            dest="overwrite",
        )
        pipeline_arg_groups["optional"].add_argument(
            "-h", "--help", action="help", help="show this help message and exit"
        )

    def returnArgs(self):
        return vars(self._pipeline_parser.parse_args())

    @classmethod
    def create(cls, config_pipelines):
        return cls(config_pipelines)

    @staticmethod
    def _yeildGroup(group_name, group_args):
        if group_name in group_args:
            for group_name, group_dict in group_args[group_name].items():
                yield group_name, group_dict

    @staticmethod
    def _processArgType(arg_args):
        try:
            arg_args["type"] = locate(arg_args["type"])
        except:
            logging.exception(f"Unable to locate type: {arg_args['type']}")
            raise

        return arg_args

    @staticmethod
    def _processArgAction(arg_args):
        # Assign the confirmDir action
        if arg_args["action"] == "confirmDir":
            arg_args["action"] = confirmDir()

        # Assign the confirmFile action
        elif arg_args["action"] == "confirmFile":
            arg_args["action"] = confirmFile()

        return arg_args

    @staticmethod
    def _processArgDefaults(arg_args):
        # Skip if the default is not a dictionary
        if not isinstance(arg_args["default"], dict):
            return arg_args

        # Assign the default string
        default_str = arg_args["default"]["str"]

        # Process suffix strings
        if "suffix" in arg_args["default"]:
            # Check if the suffix is a string, then convert to a list
            if isinstance(arg_args["default"]["suffix"], str):
                arg_args["default"]["suffix"] = [arg_args["default"]["suffix"]]

            # Loop the suffixes
            for suffix in arg_args["default"]["suffix"]:
                # Add underscores if needed
                if default_str:
                    default_str += "_"

                # Process suffix strings
                if isinstance(suffix, str):
                    default_str += suffix

                # Process suffix dicts
                elif isinstance(suffix, dict):
                    # Convert the suffix dict to lowercase keys
                    suffix_dict = {_k.lower(): _v for _k, _v in suffix.items()}

                    if "function" not in suffix_dict:
                        raise Exception(f"Suffix dict not supported: {suffix_dict}")

                    if suffix["function"] == "jobTimeStamp":
                        default_str += jobTimeStamp()
                    elif suffix["function"] == "jobRandomString":
                        default_str += jobRandomString()
                    else:
                        raise Exception(f"Function not supported: {suffix['function']}")

        # Assign the default string
        arg_args["default"] = default_str

        return arg_args

    @staticmethod
    def _processArgWildcard(pipeline_arg_name, arg_args, arg_wildcards):
        # Check if the wildcards is valid
        if arg_args["wildcards"] not in arg_wildcards:
            raise Exception(
                f'Wildcards not found ({arg_args["wildcards"]}) for argument {pipeline_arg_name}'
            )

        # Confirm no default value is present
        if "default" in arg_args:
            raise Exception(
                f'Wildcards are incompatible with default values ({arg_args["default"]}) for argument {pipeline_arg_name}'
            )

        # Add the default string to the help message
        arg_args["help"] += (
            f' (supported wildcards: {", ".join(arg_wildcards[arg_args["wildcards"]])})'
        )

        del arg_args["wildcards"]

        return arg_args


# Create variables to store randomString and timeStamp
random_string = None
time_stamp = None
