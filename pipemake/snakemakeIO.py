import os
import re
import yaml
import shutil
import logging

from pathlib import Path
from collections import defaultdict

from pipemake.singularityIO import Singularity


class SnakePipelineIO:
    def __init__(
        self,
        workflow_dir="",
        pipeline_storage_dir="",
        pipeline_job_dir="",
        singularity_dir="",
        resource_yml=None,
        scale_threads=0,
        scale_mem=0,
        indent_style="\t",
        overwrite=True,
        **kwargs,
    ):
        # Assign the basic arguments
        self._workflow_dir = workflow_dir

        # Assign the workflow filenames
        self.smkp_filename = "Snakefile"
        self.smkp_config_yml = "config.yml"
        self.smkp_resource_yml = "resources.yml"

        # Check if the overwrite is a bool or None
        if not isinstance(overwrite, bool):
            raise Exception(f"Invalid overwrite type: {overwrite}")

        # Confirm the file exists
        if (
            os.path.isfile(os.path.join(self._workflow_dir, self.smkp_filename))
            and not overwrite
        ):
            raise IOError(
                f"SnakePipeline file already exists: {os.path.join(self._workflow_dir, self.smkp_filename)}. Please rename the SnakePipeline file or overwrite"
            )

        # Confirm the resource yml is bool or None
        if not isinstance(resource_yml, (bool, type(None))):
            raise Exception(f"Invalid resource yml type: {resource_yml}")

        # Confirm the scale threads and mem are integers
        if not isinstance(scale_threads, float):
            raise Exception(f"Invalid scale threads type: {scale_threads}")
        if not isinstance(scale_mem, float):
            raise Exception(f"Invalid scale mem type: {scale_mem}")

        # Confirm the indent style
        if indent_style not in [" ", "\t"]:
            raise Exception(f"Specified indent style not supported: {indent_style}")

        # Assign the basic arguments
        self._pipe_file = open(
            os.path.join(self._workflow_dir, self.smkp_filename), "w"
        )
        self._resource_yml = resource_yml
        self._scale_threads = scale_threads
        self._scale_mem = scale_mem
        self._mem_types = [
            "mem_b",
            "mem_kb",
            "mem_mb",
            "mem_gb",
            "mem_tb",
            "mem_pb",
            "mem_kib",
            "mem_mib",
            "mem_gib",
            "mem_tib",
            "mem_pib",
        ]
        self._indent_style = indent_style

        # Assign the module storage directory, and confirm it exists
        self._module_storage_dir = os.path.join(pipeline_storage_dir, "modules")
        if not os.path.exists(self._module_storage_dir):
            raise IOError(
                f"Unable to open: {self._module_storage_dir}. Please confirm the directory exists."
            )

        # Assign the script storage directory, and confirm it exists
        self._script_storage_dir = os.path.join(pipeline_storage_dir, "scripts")
        if not os.path.exists(self._script_storage_dir):
            raise IOError(
                f"Unable to open: {self._script_storage_dir}. Please confirm the directory exists."
            )

        # Assign the module config directory
        self._internal_module_job_dir = os.path.join(pipeline_job_dir, "modules")
        self._module_job_dir = os.path.join(
            self._workflow_dir, self._internal_module_job_dir
        )
        if not os.path.exists(self._module_job_dir):
            os.makedirs(self._module_job_dir)

        # Assign the script config directory
        self._script_job_dir = os.path.join(
            self._workflow_dir, pipeline_job_dir, "scripts"
        )

        # Assign the backup directory
        self._backup_dir = os.path.join(self._workflow_dir, pipeline_job_dir, "backups")
        if not os.path.exists(self._backup_dir):
            os.makedirs(self._backup_dir)

        # Create arguments to populate
        self._module_filenames = []
        self._linked_rules = {}
        self._output_list = []
        self._config_params = set()
        self._resource_params = defaultdict(lambda: defaultdict(int))
        self._singularity_dir = singularity_dir
        self._pipeline_singularity_dict = {}
        self._pipeline_IO_attributes = {}
        self._pipeline_rule_orders = []

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    @classmethod
    def open(cls, *args, **kwargs):
        return cls(*args, **kwargs)

    def addModule(self, module_filename, **kwargs):
        # Update the module filename, if necessary
        if not module_filename.endswith(".smk"):
            module_filename = f"{module_filename}.smk"

        # Create the filepath to the stored module
        storage_module_path = os.path.join(self._module_storage_dir, module_filename)

        # Confirm the stored module exists
        if not os.path.isfile(storage_module_path):
            raise IOError(
                f"Unable to open: {storage_module_path}. Please confirm the file exists."
            )

        # Open and process the stored module file
        smk_module = SnakeFileIO.open(
            storage_module_path, singularity_dir=self._singularity_dir
        )

        # Check if the module has script files
        if len(smk_module._file_script_files) > 0:
            # Create the script job directory
            if not os.path.exists(self._script_job_dir):
                os.makedirs(self._script_job_dir)

        # Loop the rules to assign the output
        for rule in smk_module._file_rules:
            # Skip the rule if not the output rule
            if rule.rule_name != "all":
                continue

            # Check if the rule has output
            if rule._rule_output_list:
                # Save the output list from the rule
                for rule_output in rule._rule_output_list:
                    self._output_list.append(
                        _convert_indent(
                            rule_output, rule._indent_style, self._indent_style
                        )
                    )

        # Loops the file rules in test
        for rule in smk_module._file_rules:
            # Save the config params from the rule
            self._config_params.update(rule._rule_config_params)

            # Scale and save the resource params
            for resource_type, resource_value in rule._rule_resource_params.items():
                if resource_type == "threads" and self._scale_threads:
                    resource_value = self._scaler(resource_value, self._scale_threads)
                elif (
                    "mem_" in resource_type
                    and resource_type.lower() in self._mem_types
                    and self._scale_mem
                ):
                    resource_value = self._scaler(resource_value, self._scale_mem)
                self._resource_params[rule.rule_name][resource_type] = resource_value

        # Create the job module filepath
        storage_module_path = os.path.join(self._module_job_dir, module_filename)

        # Create the job module file
        smk_module.write(storage_module_path)

        # Add the module filename to the list
        self._module_filenames.append(module_filename)

        # Assign the snakefile IO attributes to the pipeline
        for rule, rule_io_attributes in smk_module._file_IO_attributes.items():
            if rule not in self._pipeline_IO_attributes:
                self._pipeline_IO_attributes[rule] = rule_io_attributes
            else:
                raise Exception(f"Duplicate rule found in pipeline: {rule}")

        # Assign the rule orders to the pipeline
        for rule_order in smk_module._rule_order_list:
            if rule_order not in self._pipeline_rule_orders:
                self._pipeline_rule_orders.append(rule_order)

        # Store the singularity containers from the module
        for (
            singularity_path,
            singularity_container,
        ) in smk_module._file_singularity_dict.items():
            if singularity_path not in self._pipeline_singularity_dict:
                self._pipeline_singularity_dict[singularity_path] = (
                    singularity_container
                )

        for script_file in smk_module._file_script_files:
            shutil.copy(
                os.path.join(self._script_storage_dir, script_file),
                os.path.join(self._script_job_dir, script_file),
            )
            logging.info(f"Script added to pipeline: {script_file}")

        logging.info(f"Module added to pipeline: {module_filename}")

    def addSnakeLink(self, link_statement):
        # Assign the input and output rules from the link statement
        input_rule = link_statement["input"]
        output_block = link_statement["output"]

        # Create list to store the linked files
        file_mappings = []

        # Assign input and output rules if specified
        if "files" in link_statement:
            for files in link_statement["file_mappings"]:
                file_mappings.append(files)

        # Convert the output block to a list for consistency
        if isinstance(output_block, str):
            output_block = [output_block]

        # Loop the output rule(s) in the output block
        for output_rule in output_block:
            # Confirm the rules are not the same
            if input_rule == output_rule:
                raise Exception(f"Input and output rules are the same: {input_rule}")

            # Check if the input rule is in the pipeline IO attributes
            if input_rule not in self._pipeline_IO_attributes:
                raise Exception(
                    f"Input rule not found in pipeline IO attributes: {input_rule}"
                )

            # Check if the output rule is in the pipeline IO attributes
            if output_rule not in self._pipeline_IO_attributes:
                raise Exception(
                    f"Output rule not found in pipeline IO attributes: {output_rule}"
                )

            # Check if the rules have linked files
            if not file_mappings:
                # Confirm the rules have the same attribute names
                if not (
                    set(self._pipeline_IO_attributes[input_rule]["output"])
                    <= set(self._pipeline_IO_attributes[output_rule]["input"])
                ):
                    raise Exception(
                        f"Input and output rules have different attribute names: {input_rule}, {output_rule}"
                    )

                # Create the linked files from the input and output rules
                file_mappings = [
                    {"input": _f, "output": _f}
                    for _f in set(self._pipeline_IO_attributes[input_rule]["output"])
                    & set(self._pipeline_IO_attributes[output_rule]["input"])
                ]

            else:
                for file_mapping in file_mappings:
                    # Check if the linked file is in the input rule
                    if (
                        file_mapping["input"]
                        not in self._pipeline_IO_attributes[input_rule]["output"]
                    ):
                        raise Exception(
                            f"Linked file not found in input rule: {file_mapping}"
                        )

                    # Check if the linked file is in the output rule
                    if (
                        file_mapping["output"]
                        not in self._pipeline_IO_attributes[output_rule]["input"]
                    ):
                        raise Exception(
                            f"Linked file not found in output rule: {file_mapping}"
                        )

            # Create the link rule
            link_rule = f"use rule {output_rule} as link_{output_rule} with:\n{self._indent_style}input:\n"

            # Add the linked attributes to the link rule
            for file_mapping in file_mappings:
                link_file = "\n".join(
                    self._pipeline_IO_attributes[input_rule]["output"][
                        file_mapping["input"]
                    ]
                )
                link_rule += f"{self._indent_style}{self._indent_style}{file_mapping['output']}={link_file}\n"

            # Check if the output rule has additional arguments
            mapped_output_files = [_f["output"] for _f in file_mappings]
            unmapped_files = dict(
                (_a, _f)
                for _a, _f in self._pipeline_IO_attributes[output_rule]["input"].items()
                if _a not in mapped_output_files
            )

            # Loop the unmapped files
            for unmapped_arg, unmapped_file in unmapped_files.items():
                unmapped_file_str = "\n".join(unmapped_file)
                link_rule += f"{self._indent_style}{self._indent_style}{unmapped_arg}={unmapped_file_str}\n"

            # Save the linked rule to the list
            self._linked_rules[output_rule] = link_rule

    def buildSingularityContainers(self):
        # Loop the singularity containers
        for singularity_container in self._pipeline_singularity_dict.values():
            # Download the singularity container
            singularity_container.download()

    def buildSingularityURLTable(self):
        # Create file in the backup directory
        url_table_file = open(
            os.path.join(self._backup_dir, "singularity_urls.tsv"), "w"
        )

        # Write the header
        url_table_file.write("Singularity Path\tURL\n")

        # Loop the singularity containers
        for singularity_container in self._pipeline_singularity_dict.values():
            # Write the singularity path and URL
            url_table_file.write(
                f"{singularity_container._image_filename}\t{singularity_container._url}\n"
            )

        # Close the file
        url_table_file.close()

    def writeConfig(self, pipeline_args):
        # Create yml dicts
        yml_config_dict = {}
        yml_resource_dict = {"resources": {}}

        # Create a list of config params to sort by length
        config_params_list = list(self._config_params)
        config_params_list.sort(key=lambda t: len(t))

        # Populate the config yml dict with the pipeline args
        for config_param in config_params_list:
            # Check if the config param is within a group
            if len(config_param) > 1:
                group, config_arg = config_param

            # Assign the config arg from the config param
            else:
                group = None
                config_arg = config_param[0]

            # Create the equivalent pipeline arg for the config arg
            pipeline_arg = config_arg.replace("-", "_")

            # Raise an error if the pipeline arg is not in the pipeline args
            if pipeline_arg not in pipeline_args:
                raise Exception(f"Pipeline arg not found: {pipeline_arg}")

            # Check if the group is not in the yml dict
            if group and group not in yml_config_dict:
                yml_config_dict[group] = {}

            # Skip the pipeline arg if it is None
            if pipeline_args[pipeline_arg] is None:
                continue

            if not group:
                yml_config_dict[config_arg] = pipeline_args[pipeline_arg]
            else:
                yml_config_dict[group][config_arg] = pipeline_args[pipeline_arg]

        # Populate the resource yml dict
        for rule_name, rule_resource_dict in self._resource_params.items():
            if not rule_resource_dict:
                continue

            # Create the rule resource dict
            yml_resource_dict["resources"][rule_name] = {}

            # Populate the rule resource dict
            for resource_name, resource_value in rule_resource_dict.items():
                yml_resource_dict["resources"][rule_name][resource_name] = (
                    resource_value
                )

        # Check if the resource yml should be a separate file
        if not self._resource_yml:
            yml_config_dict.update(yml_resource_dict)
        else:
            self._createYml(
                yml_resource_dict,
                os.path.join(self._workflow_dir, self.smkp_resource_yml),
            )
            logging.info("Resource file written successfully")

        self._createYml(
            yml_config_dict, os.path.join(self._workflow_dir, self.smkp_config_yml)
        )
        logging.info("Config file written successfully")

    def writePipeline(self):
        # Assign the config file(s) and working directory
        self._pipe_file.write(f"configfile: '{self.smkp_config_yml}'\n")
        if self._resource_yml:
            self._pipe_file.write(f"configfile: '{self.smkp_resource_yml}'\n")

        # Create the rule all block
        self._pipe_file.write("\nrule all:\n")
        self._pipe_file.write(f"{self._indent_style}input:\n")

        # Add the rule all input
        indented_input = "\n".join(_output for _output in self._output_list)
        self._pipe_file.write(indented_input + "\n\n")

        # Create the input block of modules
        for module_filename in self._module_filenames:
            self._pipe_file.write(
                f'include: "{os.path.join(self._internal_module_job_dir, module_filename)}"\n'
            )

        logging.info(f"Module files written successfully: {self._module_filenames}")

        # Skip rule linking if no linked rules
        if not self._linked_rules.values():
            logging.info("Pipeline written successfully")
            return

        # Assign linked rule orders
        linked_rule_orders = self._orderLinkedRules(
            self._pipeline_rule_orders, self._linked_rules.keys()
        )

        if linked_rule_orders:
            self._pipe_file.write("\n")

        # Assign the rule orders for the linked rules
        for linked_order in linked_rule_orders:
            self._pipe_file.write(linked_order)

        self._pipe_file.write("\n")

        # Create the linked rules
        for linked_rule in self._linked_rules.values():
            self._pipe_file.write(linked_rule + "\n")

        logging.info("Pipeline written successfully")

    def close(self):
        self._pipe_file.close()

        # Create the singularity url table, if necessary
        if self._pipeline_singularity_dict:
            self.buildSingularityURLTable()

        # Create backups of the pipeline files
        shutil.copy(
            os.path.join(self._workflow_dir, self.smkp_filename),
            os.path.join(self._backup_dir, f"{self.smkp_filename}.bkp"),
        )
        shutil.copy(
            os.path.join(self._workflow_dir, self.smkp_config_yml),
            os.path.join(self._backup_dir, f"{self.smkp_config_yml}.bkp"),
        )
        if os.path.isfile(os.path.join(self._workflow_dir, self.smkp_resource_yml)):
            shutil.copy(
                f"{os.path.join(self._workflow_dir, self.smkp_resource_yml)}",
                os.path.join(self._backup_dir, f"{self.smkp_resource_yml}.bkp"),
            )

        logging.info("Pipeline backups created successfully")

    @staticmethod
    def _orderLinkedRules(rule_orders, rules_to_link):
        """
        Create a linked rule order string from the given rule orders and rules to link.
        """

        # Create list to store the linked rule orders
        link_rule_order_list = []

        # Create a dictionary to store the rule orders that contain rules being linked
        rule_orders_w_links = defaultdict(list)

        # Assign the order index of the rule string in the rule_orders
        for index, rule_order in enumerate(rule_orders):
            for rule in rules_to_link:
                if rule in rule_order:
                    rule_orders_w_links[index].append(rule)

        # Remove rule orders with fewer than 2 rules
        for rule_order_index in list(rule_orders_w_links):
            if len(rule_orders_w_links[rule_order_index]) < 2:
                del rule_orders_w_links[rule_order_index]

        # Loop through the rule orders to create the linked rule orders
        for rule_order_index, rules in rule_orders_w_links.items():
            # Assign the position of each rule in the rule_order str
            rule_order_str = rule_orders[rule_order_index]

            # Check if the order symbols are ambiguous
            if ">" in rule_order_str and "<" in rule_order_str:
                raise ValueError(
                    "Order symbols are ambiguous. Please use only one order symbol in the rule order string."
                )

            # Check if the order symbols are missing
            if ">" not in rule_order_str and "<" not in rule_order_str:
                raise ValueError(
                    "Order symbols are missing. Please use either '>' or '<' in the rule order string."
                )

            # Assign the order symbol based on the rule order string
            if ">" in rule_order_str:
                order_symbol = ">"
            elif "<" in rule_order_str:
                order_symbol = "<"
            else:
                raise Exception(f"Unexpected error: {rule_order_str}")

            # Create a string to store the linked rule order
            link_rule_order_str = "ruleorder: "

            # Create a sorted list of (rule_position, rule)
            link_rule_order = [
                (re.search(r"\b" + _r + r"\b", rule_order_str).start(), _r)
                for _r in rules
            ]
            link_rule_order.sort(key=lambda x: x[0])

            # Create the linked rule order string
            for link_rule_pos, link_rule_group in enumerate(link_rule_order):
                # Assign the first rule in the order
                if link_rule_pos == 0:
                    link_rule_order_str += f"link_{link_rule_group[1]}"
                    continue

                # Add the rule to the rule order string
                link_rule_order_str += f" {order_symbol} link_{link_rule_group[1]}"

            # Save the linked rule order string to the list
            link_rule_order_list.append(link_rule_order_str + "\n")

        return link_rule_order_list

    @staticmethod
    def _scaler(value, scaler, value_type=float, output_type=int, **kwargs):
        # Confirm a scaler was given
        if not scaler:
            raise Exception(f"No scaler given for: {value}")

        # Convert the value to the correct type
        try:
            value = value_type(value)
        except:
            logging.exception(f"Unable to convert value to type: {value}")
            raise

        # Scale and return the input
        return output_type(value * scaler)

    @staticmethod
    def _createYml(yml_dict, yml_filename, **kwargs):
        class MyDumper(yaml.Dumper):
            def increase_indent(self, flow=False, indentless=False):
                return super(MyDumper, self).increase_indent(flow, False)

        yml_file = open(yml_filename, "w")
        yml_file.write(yaml.dump(yml_dict, Dumper=MyDumper, sort_keys=False))
        yml_file.close()


class SnakeFileIO:
    def __init__(self, smk_filename, singularity_dir="", **kwargs):
        # Confirm the file exists
        if not os.path.isfile(smk_filename):
            raise IOError(
                f"Unable to open: {smk_filename}. Please confirm the file exists."
            )

        # Assign the basic arguments
        self.filename = smk_filename
        self._singularity_dir = singularity_dir
        self._indent_style = None
        self._output_rule = "all"
        self._exclude_rules = [self._output_rule]
        self._file_rules = []
        self._non_rule_text = ""

        # Assign the indent style
        self._assignIndent()

        # Parse the snakefile
        self._parseFile()

    @property
    def _rule_order_list(self):
        rule_order = []
        for text in self._non_rule_text.splitlines():
            if text.startswith("ruleorder:"):
                rule_order.append(text.split("ruleorder:")[1].strip())
        return rule_order

    @property
    def _file_IO_attributes(self):
        file_io_attributes = {}
        for rule in self._file_rules:
            if rule.rule_name == self._output_rule:
                continue

            file_io_attributes[rule.rule_name] = rule._rule_IO_attributes
        return file_io_attributes

    @property
    def _file_singularity_dict(self):
        file_singularity_dict = {}
        for rule in self._file_rules:
            for (
                singularity_path,
                singularity_container,
            ) in rule._rule_singularity_dict.items():
                if singularity_path not in file_singularity_dict:
                    file_singularity_dict[singularity_path] = singularity_container
        return file_singularity_dict

    @property
    def _file_script_files(self):
        file_script_files = []
        for rule in self._file_rules:
            for script_file in rule._rule_script_files:
                file_script_files.append(script_file)
        return file_script_files

    def _assignIndent(self):
        # Bool to store if within a rule
        in_rule_block = False

        with open(self.filename) as smk_file:
            for smk_line in smk_file:
                # Skip blank lines
                if not smk_line.strip():
                    continue

                # Confirm indent style once assigned
                if self._indent_style:
                    # Skip rules as they have no indent
                    if (
                        smk_line.startswith("rule")
                        or smk_line.startswith("checkpoint")
                        or smk_line.startswith("def")
                    ):
                        continue

                    if not smk_line.startswith(self._indent_style):
                        raise Exception(
                            f'Inconsistent indent style "{smk_line[0]}" in: {self.filename}'
                        )

                    # Get the count of the indent style
                    indent_count = len(smk_line) - len(smk_line.lstrip())

                    # Calc the indent remainder
                    indent_remainder = indent_count % len(self._indent_style)

                    # Check indent remainder
                    if indent_remainder:
                        raise Exception(f"Indent style remainder: {indent_remainder}")

                # Assign the indent style
                elif in_rule_block:
                    # Get the count of the leading whiespace
                    leading_whitespace_count = len(smk_line) - len(smk_line.lstrip())

                    # Check for table indent style
                    if smk_line[0] == "\t":
                        self._indent_style = "\t" * leading_whitespace_count
                        logging.info(
                            f"Indent style: Tab(s). Count: {leading_whitespace_count}"
                        )

                    # Check for space indent style
                    elif smk_line[0] == " ":
                        self._indent_style = " " * leading_whitespace_count
                        logging.info(
                            f"Indent style: Spaces(s). Count: {leading_whitespace_count}"
                        )

                    # Report error, if other
                    else:
                        raise Exception("Unable to assign indent style")

                # Check if within first rule block
                if smk_line.split(" ")[0] == "rule":
                    in_rule_block = True

    def _parseFile(self):
        # Create a string to store the rule block
        rule_block = ""

        # Create bool to check if in the config block
        within_rule_block = False

        # Parse the snakefile
        with open(self.filename) as smk_file:
            for smk_line in smk_file:
                # Skip blank lines
                if not smk_line.strip():
                    continue

                # Split the line by the indent style (rule, attribute, param)
                smk_list = smk_line.rstrip().split(self._indent_style)

                # Assign if line is a rule or checkpoint
                is_rule_line = (
                    smk_list[0].startswith("rule")
                    or smk_list[0].startswith("checkpoint")
                ) and smk_line.rstrip().endswith(":")

                # Check if the line is the end of a rule
                if smk_list[0] and not is_rule_line:
                    within_rule_block = False

                # Check if the line is the start of a rule
                if is_rule_line:
                    within_rule_block = True

                    # Check if a previous rule block was found
                    if rule_block:
                        self._file_rules.append(
                            SnakeRuleIO.read(
                                rule_block,
                                singularity_dir=self._singularity_dir,
                                indent_style=self._indent_style,
                                output_rule=self._output_rule,
                            )
                        )
                        rule_block = ""

                # Add the line to the rule block
                if within_rule_block and smk_line.strip():
                    rule_block += smk_line

                # Save text that is not within a rule
                if not within_rule_block and smk_line.strip():
                    self._non_rule_text += smk_line

            # Parse the final rule block
            if rule_block:
                self._file_rules.append(
                    SnakeRuleIO.read(
                        rule_block,
                        singularity_dir=self._singularity_dir,
                        indent_style=self._indent_style,
                        output_rule=self._output_rule,
                    )
                )

    def write(self, filename):
        # Open the file
        with open(filename, "w") as smk_file:
            # Check if there is any non-rule text
            if self._non_rule_text:
                smk_file.write(self._non_rule_text)

            # Loop the file rules
            for rule in self._file_rules:
                # Write the rule, if included in the output
                if not rule._output_rule:
                    smk_file.write(f"\n\n{rule}")

    @classmethod
    def open(cls, *args, **kwargs):
        return cls(*args, **kwargs)


class SnakeRuleIO:
    def __init__(
        self, rule_str, singularity_dir="", indent_style=None, output_rule="", **kwargs
    ):
        # Assign the rule argument to populate
        self._original_text = rule_str
        self._original_text_list = self._original_text.splitlines()
        self.rule_name = self._original_text_list[0].split()[1][:-1]
        self._output_rule = self.rule_name == output_rule
        self._singularity_dir = singularity_dir
        self._indent_style = indent_style
        self._rule_config_params = set()
        self._rule_singularity_dict = {}
        self._rule_script_files = []
        self._rule_resource_params = {}
        self._rule_text = self._original_text_list[0] + "\n"
        self._rule_output_list = []
        self._rule_IO_attributes = defaultdict(
            list, {k: defaultdict(list) for k in ("input", "output")}
        )

        # Assign the fixed arguments
        self._resource_attributes = ["threads", "resources"]

        # Update and assign the rule parameters
        if self._output_rule:
            self._setOutput()
            self._setParams()

        # Set the output list
        elif not self._output_rule:
            self._updateRule()
            self._setParams()
            self._assignIOAttributes()

    def __str__(self):
        return self._rule_text

    def __repr__(self):
        return self._rule_text

    @classmethod
    def read(cls, rule_str, *args, **kwargs):
        return cls(rule_str, *args, **kwargs)

    def _updateResource(self, resource_attribute, rule_line, attribute_level):
        # Assign the resource str
        resource_str = rule_line.strip()

        # Check if the resource is a plain integer
        if resource_str.isdigit():
            resource_name = resource_attribute
            resource_value = int(resource_str)

        else:
            # Assign the expected resource syntax
            resource_syntax_dict = {1: ":", 2: "="}

            # Check if the resource attributes can be assigned
            if resource_syntax_dict[attribute_level] not in rule_line:
                raise Exception(
                    f"Resource attribute error for {resource_attribute}. Value not assigned: {rule_line}"
                )

            # Assign the resource name and value
            resource_name, resource_value = rule_line.strip().split(
                resource_syntax_dict[attribute_level]
            )

            # Remove any whitespaces
            resource_name = resource_name.strip()
            resource_value = resource_value.strip()

            # Remove any trailing commas
            if resource_value.endswith(","):
                resource_value = resource_value[:-1]

            # Check if the resource defines the shell executable, and return the rule line if so
            if resource_name == "shell_exec":
                self._rule_resource_params[resource_name] = resource_value.replace(
                    '"', ""
                ).replace("'", "")
                return rule_line

            # Determine if the resource value type
            if not resource_value.isdigit():
                raise Exception(
                    f"Complex resource support not yet implemented: {rule_line}"
                )

        # Assign the resource value
        self._rule_resource_params[resource_name] = int(resource_value)

        # Assign the resource config attribute
        resource_config = f'config["resources"]["{self.rule_name}"]["{resource_name}"]'

        # Update the resource list
        rule_line = rule_line.replace(str(resource_value), resource_config)

        return rule_line

    def _updateContainer(self, container_str):
        # Check if a singularity directory is specified
        if not self._singularity_dir:
            return container_str

        # Assign the container text
        container_url = container_str.strip().split('"')[1]

        # Store the container
        container = Singularity.fromURL(container_url, self._singularity_dir)

        # Add the container to the singularity dict
        self._rule_singularity_dict[container.returnPath()] = container

        return container.updateContainer(container_str)

    def _updateScript(self, script_str, script_dir="scripts"):
        # Assign the script path
        script_path = Path(script_str.strip().split('"')[1])

        # Add the script file to the list
        self._rule_script_files.append(script_path.name)

        # Assign the script path
        script_output_path = os.path.join("..", script_dir, script_path.name)

        # Return the container path
        return script_str.split('"')[0] + f'"{script_output_path}"'

    def _setParams(self):
        def mergeLines(str):
            return "".join([_s.strip() for _s in str.splitlines()])

        # Find all configs, separated by optional and required
        for config_match in re.finditer(
            r"(?<=config)(\[[^:\]]+\]+)+|if\s*(.*?)\s*in\s*?config((\[(?![^\]]*:[^\]]*)[^\]]*\])*)?",
            self._original_text,
        ):
            # Check if the match is an optional config param
            if "if " in config_match.group(0):
                config_str = f"[{config_match.group(2)}]"
                if config_match.group(3):
                    config_str = mergeLines(config_match.group(3)) + config_str

            # Check if the match is a required config param
            else:
                config_str = mergeLines(config_match.group(0))

            # Remove quotes from the config string
            config_str = config_str.replace('"', "").replace("'", "")

            # Split the param str by either square brackets using regex
            config_list = tuple([_p for _p in re.split(r"\[|\]", config_str) if _p])

            # Check if the match is associated with a value
            config_has_value = (
                False if self._original_text.find(f"{config_str}:") == -1 else True
            )

            # Check if the param does not have a value, i.e. config argument
            if config_has_value:
                raise Exception(f"Config param error. Has value: {self._original_text}")

            # Confirm the config list isn't more than two parameters
            if len(config_list) > 2:
                raise Exception(
                    f"Config param error. Too many param levels: {self._original_text}"
                )

            # Add the config assignment to the set
            self._rule_config_params.add(tuple(config_list))

    def _updateRule(self):
        for rule_line in self._original_text_list[1:]:
            if not rule_line.strip():
                continue

            # Split the line by the indent style
            rule_line_list = rule_line.rstrip().split(self._indent_style)

            # Assign the attribute level
            attribute_level = _attributeLevel(rule_line_list)

            # Check if the line is a rule attribute
            if attribute_level == 1:
                attribute_name, attribute_value = [
                    _a.strip() for _a in rule_line_list[1].split(":", 1)
                ]

            # Process resource attributes
            if attribute_name in self._resource_attributes:
                if attribute_value:
                    rule_line = self._updateResource(
                        attribute_name, rule_line, attribute_level
                    )
                elif attribute_level > 1:
                    rule_line = self._updateResource(
                        attribute_name, rule_line, attribute_level
                    )

            # Process container attributes
            elif attribute_name == "singularity":
                if attribute_value:
                    rule_line = self._updateContainer(rule_line)
                elif attribute_level > 1:
                    rule_line = self._updateContainer(rule_line)

            # Process script attributes
            elif attribute_name == "script":
                if attribute_value:
                    rule_line = self._updateScript(rule_line)
                elif attribute_level > 1:
                    rule_line = self._updateScript(rule_line)

            # Update the rule text
            self._rule_text += rule_line + "\n"

    def _setOutput(self):
        for rule_line in self._original_text_list[1:]:
            if not rule_line.strip():
                continue

            # Split the line by the indent style
            rule_line_list = rule_line.rstrip().split(self._indent_style)

            # Assign the attribute level
            attribute_level = _attributeLevel(rule_line_list)

            # Check if the line is a rule attribute
            if attribute_level == 1 and rule_line.rstrip().endswith(":"):
                attribute_name = rule_line_list[1][:-1]

            # Check if the rule is output, confirm the attribute is input, and update the rule output
            elif attribute_name == "input":
                self._rule_output_list.append(rule_line.rstrip())

        # Check that output list was populated
        if not self._rule_output_list:
            raise Exception("Output list not populated for rule")

    def _assignIOAttributes(self):
        def _removeTempStatement(link_statement):
            if len(link_statement) == 1 and isinstance(link_statement[0], int):
                return link_statement
            elif "temp(" not in link_statement[0]:
                return link_statement
            elif len(link_statement) == 1:
                temp_end_pos = link_statement[0].rfind(")")
                link_statement[0] = (
                    link_statement[0][:temp_end_pos]
                    + link_statement[0][temp_end_pos + 1 :]
                )
                link_statement[0] = link_statement[0].replace("temp(", "")
                return link_statement
            elif len(link_statement) > 1:
                link_statement = link_statement[1:-1]
                if not link_statement[-1].endswith(","):
                    link_statement[-1] = link_statement[-1] + ","
                link_statement = [
                    _s[len(self._indent_style) :] for _s in link_statement
                ]
                link_statement[0] = link_statement[0].strip()
                return link_statement

        # Create string to store the attribute id, as it be on multiple lines
        attribute_name = ""

        # Loop through the rule list
        for rule_line in self._original_text_list[1:]:
            if not rule_line.strip():
                continue

            # Split the line by the indent style
            rule_line_list = rule_line.rstrip().split(self._indent_style)

            # Assign the attribute level
            attribute_level = _attributeLevel(rule_line_list)

            # Check if the line is a rule attribute
            if attribute_level == 1 and rule_line.rstrip().endswith(":"):
                attribute_group = rule_line_list[1][:-1]

            # Check if the rule is output, confirm the attribute is input, and update the rule output
            elif attribute_group in ["input", "output"]:
                if attribute_level == 2:
                    # Check if any alphanumeric characters, if not just append the attribute text (multi-line attributes)
                    if not any(_c.isalnum() for _c in rule_line_list[2]):
                        self._rule_IO_attributes[attribute_group][
                            attribute_name
                        ].append(
                            _attributeStr(
                                rule_line_list[2], attribute_level, self._indent_style
                            )
                        )

                    # Assign the attribute with a known name
                    elif "=" in rule_line_list[2]:
                        attribute_name, attribute_str = rule_line_list[2].split("=")
                        self._rule_IO_attributes[attribute_group][
                            attribute_name
                        ].append(attribute_str)

                    # Assign the attribute by order, as name is not known
                    else:
                        attribute_name = len(self._rule_IO_attributes[attribute_group])
                        attribute_str = rule_line_list[2]

                        self._rule_IO_attributes[attribute_group][
                            attribute_name
                        ].append(attribute_str)

                # Append attribute text to the current attribute name (multi-line attributes)
                elif attribute_level > 2:
                    self._rule_IO_attributes[attribute_group][attribute_name].append(
                        _attributeStr(
                            rule_line_list[attribute_level],
                            attribute_level,
                            self._indent_style,
                        )
                    )

        for io_group, group_attribute_dict in self._rule_IO_attributes.items():
            for group_attribute, group_value in group_attribute_dict.items():
                self._rule_IO_attributes[io_group][group_attribute] = (
                    _removeTempStatement(group_value)
                )


def _attributeLevel(attribute_list):
    for attribute_level, attribute_item in enumerate(attribute_list):
        if not attribute_item:
            continue
        return attribute_level


def _attributeStr(attribute_str, attribute_level, indent_style):
    # Check if the attribute level is valid
    if attribute_level < 0:
        raise ValueError("Invalid attribute level. Must be greater than or equal to 0.")

    return (indent_style * attribute_level) + attribute_str


def _convert_indent(source_str, source_indent_style, target_indent_style):
    # Get the count of the source indent style
    indent_count = source_str.count(source_indent_style)
    # Return the converted indent style
    return f"{indent_count*target_indent_style}{source_str.strip()}"
