import os
import re
import yaml
import logging

from pipemake.processIO import processInput


class ConfigPipelinesIO(dict):
    def __init__(self, *arg, **kwargs):
        super(ConfigPipelinesIO, self).__init__(*arg, **kwargs)

    @classmethod
    def fromDirectory(cls, directory):
        # Check the config directory exists
        config_dir = os.path.join(directory, "configs")
        if not os.path.isdir(config_dir):
            raise Exception(f"Unable to find pipeline directory: {config_dir}")

        # Create dict to store the pipeline configs
        pipeline_configs = {}

        # Loop the configs
        for config_filename in os.listdir(config_dir):
            config_file = os.path.join(config_dir, config_filename)

            # Confirm the path is a file
            if not os.path.isfile(config_file):
                continue

            # Load the config file
            pipeline_config = ConfigPipelineIO.fromYAML(config_file)

            # Store the pipeline config
            pipeline_configs[pipeline_config.name] = pipeline_config

        return cls(**pipeline_configs)


class ConfigPipelineIO:
    def __init__(self, config_dict):
        # Check if the pipeline name is provided
        if "pipeline" not in config_dict:
            raise Exception("No pipeline name provided in the configuration file")

        # Set the pipeline name
        self.name = config_dict["pipeline"]
        config_dict.pop("pipeline")

        # Check if the pipeline version is provided
        if "version" not in config_dict:
            raise Exception("No pipeline version provided in the configuration file")

        # Set the pipeline version
        self.version = config_dict["version"]
        config_dict.pop("version")

        # Check if the pipeline parser is provided
        if "parser" not in config_dict:
            raise Exception("No parser provided in the configuration file")

        # Set the pipeline parser
        self.parser_dict = config_dict["parser"]
        config_dict.pop("parser")

        # Check if the pipeline setup is provided
        if "setup" not in config_dict:
            logging.warning("No setup provided in the configuration file")

        # Set the pipeline setup
        self._setup_dict = config_dict["setup"]
        config_dict.pop("setup")

        # Check if the pipeline snakefiles are provided
        if "snakemake" not in config_dict:
            raise Exception("No snakemake section provided in the configuration file")

        # Check if the pipeline snakefiles are provided
        if "modules" not in config_dict["snakemake"]:
            raise Exception("No snakefiles provided in the configuration file")

        # Set the pipeline snakefiles
        self._snakefiles = set(config_dict["snakemake"]["modules"])
        logging.info(f"Loaded the following snakefiles: {self._snakefiles}")

        # Check if linked rules are provided
        if "links" in config_dict["snakemake"]:
            self._snakelinks = set(config_dict["snakemake"]["links"])
        else:
            self._snakelinks = set()

        config_dict.pop("snakemake")

        # Record unused configuration parameters
        for key in config_dict:
            logging.info(f"Unused configuration parameter: {key}")

        # Assign the other configuration arguments
        self._singularity_bindings = set()
        self._inputlinks = set()
        self.samples = {}

    @property
    def snakefiles(self):
        return list(self._snakefiles)

    @property
    def snakelinks(self):
        return list(self._snakelinks)

    @property
    def singularity_bindings(self):
        return list(self._singularity_bindings)

    @classmethod
    def fromYAML(cls, config_filename):
        with open(config_filename, "r") as config_file:
            config_dict = yaml.safe_load(config_file)
        return cls(config_dict)

    def helpMessage(self, pipeline_args):
        if not pipeline_args["singularity_dir"]:
            logging.warn(
                "No singularity path specified, this may result in redundant singularity containers. Please consider specifying the singularity path by including the --singularity-dir argument or setting the PM_SINGULARITY_DIR environment variable."
            )

        print(
            f"{self.name} version {self.version} has been configured, please use the following command within the {pipeline_args['workflow_dir']} directory:"
        )
        logging.info(
            f"{self.name} version {self.version} has been configured, please use the following command within the {pipeline_args['workflow_dir']} directory:"
        )

        if self.singularity_bindings:
            print(
                f"snakemake --use-singularity --singularity-args '--bind {','.join(self.singularity_bindings)}'"
            )
            logging.info(
                f"snakemake --use-singularity --singularity-args '--bind {','.join(self.singularity_bindings)}'"
            )
        else:
            print("snakemake --use-singularity")
            logging.info("snakemake --use-singularity")

    def setupPipeline(self, pipeline_args):
        def stringToArgs(arg_string):
            return re.findall(r"\{(.*?)\}", arg_string.replace("-", "_"))

        def processSetupArgs(arg_string):
            for setup_value_arg in stringToArgs(arg_string):
                # Confirm the setup value argument is specified
                if setup_value_arg not in pipeline_args:
                    raise Exception(
                        f"Setup argument {setup_value_arg} not found among pipeline arguments"
                    )

            # Replace the setup value arguments with the pipeline arguments
            return arg_string.replace("-", "_").format(**pipeline_args)

        # Confirm the workflow prefix is specified
        if "workflow_dir" not in pipeline_args:
            raise Exception("Workflow prefix not found among pipeline arguments")

        for setup_name, setup_dict in self._setup_dict.items():
            # Create a dict to store the setup arguments
            setup_args = {
                "method": "",
                "args": {"workflow_dir": pipeline_args["workflow_dir"]},
            }

            # Confirm the setup method was assigned
            for setup_method, assignment_arg_str in setup_dict["methods"].items():
                assignment_arg = stringToArgs(assignment_arg_str)

                # Confirm a single assignment argument is specified
                if len(assignment_arg) != 1:
                    raise Exception(
                        f"Setup method {setup_method} has more than one assignment argument, please check the configuration file"
                    )
                assignment_arg = assignment_arg[0]

                # Confirm the assignment argument is specified
                if assignment_arg not in pipeline_args:
                    raise Exception(
                        f"Setup method {assignment_arg} not found among pipeline arguments"
                    )

                # Confirm a single assignment argument was specified
                if pipeline_args[assignment_arg] is not None:
                    if setup_args["method"]:
                        raise Exception(
                            f"Setup ({setup_name}) already assigned to {setup_args['method']}, cannot assign {setup_method}"
                        )
                    else:
                        setup_args["method"] = setup_method
                        setup_args["args"][setup_method] = processSetupArgs(
                            assignment_arg_str
                        )

            # Standardize the input arguments
            for setup_arg, setup_value in setup_dict["args"].items():
                if isinstance(setup_value, str):
                    setup_value = processSetupArgs(setup_value)

                elif isinstance(setup_value, list):
                    setup_value = [processSetupArgs(_arg) for _arg in setup_value]

                # Assign the setup arguments
                setup_args["args"][setup_arg] = setup_value

            # Process the input
            processed_input = processInput(**setup_args)

            # Standardize the input
            processed_input.standardize()

            # Assign the input paths
            input_paths = processed_input.returnPaths()

            # Check for method paths, and update the singularity bindings if found
            if len(input_paths) > 0:
                self._singularity_bindings.update(input_paths)

            # Check if any sample keywords were provided
            if "sample_keywords" in setup_dict["args"]:
                input_samples = processed_input.returnSamples()

                # Confirm the samples are not already assigned
                if len(self.samples) > 0:
                    raise Exception("Samples already assigned")

                # Store the samples
                self.samples = input_samples

            if "snakefiles" in setup_dict:
                self._snakefiles.update(setup_dict["snakefiles"])

            """

            

            method_samples = returnSamples(**samples_args)

                    # Confirm the samples are not already assigned
                    if method_samples and len(self.samples) > 0:
                        raise Exception("Samples already assigned")

                    # Store the samples
                    self.samples = method_samples

                if "snakefiles" in method_args:
                    # Add method snakefiles to the pipeline
                    self._snakefiles.update(method_args["snakefiles"])

            


           
            # Create a string to store the assigned method
            assigned_method = ""

            for method_name, method_args in setup_methods.items():
                # Assign the input args
                input_args = method_args["input"]

                # Check for missing arguments
                for input_arg in input_args["args"]:
                    if input_arg.replace("-", "_") not in pipeline_args:
                        raise Exception(
                            f"Setup argument {input_arg} not found among pipeline argument"
                        )

                # Confirm expected args were specified
                method_missing_args = [
                    _a
                    for _a in input_args["args"]
                    if not pipeline_args[_a.replace("-", "_")]
                ]

                # Skip if missing arguments
                if method_missing_args:
                    continue

                # Check if the method was already assigned and report an error if so
                if assigned_method:
                    raise Exception(
                        f"Setup {setup_name} already assigned to {assigned_method}, cannot assign {method_name}"
                    )
                else:
                    assigned_method = method_name

                if "standardize" in method_args:
                    # Create a dict of the standardize args
                    std_args = copy.deepcopy(method_args["standardize"])

                    # Add the workflow prefix to the args
                    std_args["args"]["workflow_dir"] = "{workflow_dir}"

                    # Update the arguments
                    for std_arg, arg_params in std_args["args"].items():
                        if not isinstance(arg_params, str):
                            continue
                        std_args["args"][std_arg] = arg_params.replace("-", "_").format(
                            **pipeline_args
                        )

                    # Standardize the input
                    standardizeInput(**std_args)

                    # Assign the method paths
                    method_paths = returnPaths(**std_args)

                    # Check for method paths, and update the singularity bindings if found
                    if len(method_paths) > 0:
                        self._singularity_bindings.update(method_paths)

                if "samples" in method_args:
                    # Create a dict of the standardize args
                    samples_args = copy.deepcopy(method_args["samples"])

                    # Update the arguments
                    for samples_arg, arg_params in samples_args["args"].items():
                        if not isinstance(arg_params, str):
                            continue
                        samples_args["args"][samples_arg] = arg_params.replace(
                            "-", "_"
                        ).format(**pipeline_args)

                    # Assign the samples from the method
                    method_samples = returnSamples(**samples_args)

                    # Confirm the samples are not already assigned
                    if method_samples and len(self.samples) > 0:
                        raise Exception("Samples already assigned")

                    # Store the samples
                    self.samples = method_samples

                if "snakefiles" in method_args:
                    # Add method snakefiles to the pipeline
                    self._snakefiles.update(method_args["snakefiles"])
            """
