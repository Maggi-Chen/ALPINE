#!/usr/bin/env python3
"""
Skeleton of a generic launcher that can be applied to most pipelines
"""
import argparse
import logging
import sys
import os
import time
import re

# Logging Config
from log_config import enable_logging

# Pipeline Config
import pipeline_config

from samplesheet import get_samples, generate_output_file_list

# CWL Platform
from cwl_platform import SUPPORTED_PLATFORMS, PlatformFactory

# Set up logging
enable_logging(log_file=os.environ.get("LOG_FILE", "launcher.log"))
logger = logging.getLogger(__name__)


def copy_workflows(platform_config, platform, project):
    """
    Copy all workflows as specific in the pipeline config for the platform to the destination project
    :return: Dictionary of workflow -> id
    """
    workflows = {}
    for workflow in platform_config:
        if workflow.endswith("_workflow"):
            workflows[workflow] = platform.copy_workflow(
                platform_config[workflow], project
            )
    return workflows


def get_default_parameters(args):
    ''' Set default parameters '''
    parameters = {
        "reffile": {
            "class": "File",
            "path": args.reffile_id,
        },
        "primer_1": args.primer_1,
        "primer_2": args.primer_2,
        "config": {
            "class": "File",
            "path": args.config_id,
        },
        "truncated_cutoff": args.truncated_cutoff,
        "min_qual": args.min_qual,
        "five_prime_HA_seq": args.five_prime_HA_seq,
        "three_prime_HA_seq": args.three_prime_HA_seq,
        "HA_Seq_match_ratio": args.HA_match_ratio,
    }
    if args.left_itr_seq:
        parameters["left_itr_seq"] = args.left_itr_seq
    if args.right_itr_seq:
        parameters["right_itr_seq"] = args.right_itr_seq
    if args.seed_size:
        parameters["seed_size"] = args.seed_size
    if args.perc_identity:
        parameters["perc_identity"] = args.perc_identity
    return parameters


def rename_output(platform, project, merge_workflow):
    """
    Rename output files from merge steps to keep latest results without
    _[0-9]_ prefix.
    """
    for output in platform.get_task_outputs(merge_workflow):
        try:
            output_filename = platform.get_task_output_filename(merge_workflow, output)
            rename_output_file(platform, project, output_filename)
        except:
            logger.error("Failed to rename file %s", output_filename)


def rename_output_file(platform, project, filename):
    """
    Rename a file in SBG if filename starts with _[0-9]_ by
    switching filename between latest output file and previous
    output file without _[0-9]_ prefix.
    """
    if not re.match("^_[0-9][0-9]*_", filename):
        return True
    fileid = platform.get_file_id(project, filename)
    if not fileid:
        logger.warning("Cannot find file %s", filename)
        return True
    target_name = "_".join(filename.split("_")[2:])
    oldfile = platform.get_file_id(project, target_name)
    platform.rename_file(fileid, "temporary_file")
    platform.rename_file(oldfile, filename)
    platform.rename_file(fileid, target_name)
    return True


def run_merge_count_table_workflow(samples, workflow, platform, project):
    """ Run merge step to merge per-sample count tables """
    parameters = {
        "count_tables": [],
    }

    for _, sample in enumerate(samples):
        if samples[sample]["workflows"]["per_sample"]["state"] == "Complete":
            per_sample_output = platform.get_task_output(
                samples[sample]["workflows"]["per_sample"]["workflow"],
                "count_table",
            )
            if per_sample_output:
                parameters["count_tables"] += [
                    {
                        "class": "File",
                        "path": per_sample_output,
                    }
                ]
    merge_workflow = platform.submit_task("Merge count table", project, workflow, parameters)
    if merge_workflow is None:
        logger.error("Failed to submit task for Merge count table")
        return None
    return merge_workflow


def run_persample_workflow(
    samples, workflow, parameters, workflow_name, platform, project, no_reuse
):
    """
    Submits per sample workflow on the platform for all samples in the samples dictionary.

    :param samples: Dictionary of samples where sample name is dictionary key
    :param workflow: Workflow to run
    :param parameters: Parameters for the worflow
    :param workflow_name: Name of the workflow
    :param platform: Platform to run on
    :param project: Project to run in
    :return: Dictionary of samples with platform workflow info added.

    This function really is to set
        samples[sample]['state'] = 'Submitted'
        samples[sample]['workflow'] = sample_workflow

    To accomodate downstream workflows, we should expand this to
        samples[sample]['workflows']['workflow_name']['state'] = 'Submitted'
        samples[sample]['workflows']['workflow_name']['workflow'] = sample_workflow

    """
    # submit per-sample workflow for all samples on the platform
    for idx, sample in enumerate(samples):
        logger.info(
            "[%d/%d] Processing %s for %s", idx + 1, len(samples), workflow_name, sample
        )

        # If the sample doesn't have a workflows key, add it.
        if "workflows" not in samples[sample]:
            samples[sample]["workflows"] = {workflow_name: {}}
        else:
            samples[sample]["workflows"][workflow_name] = {}

        parameters["fastqfile"] = {
            "class": "File",
            "path": platform.get_file_id(
                project,
                samples[sample]["fastq"],
            ),
        }
        parameters["sample_name"] = sample

        taskname = sample + " - Per-sample Transgene analysis"
        # Get any previous workflow executions for this sample
        if no_reuse:
            tasks = None
        else:
            tasks = platform.get_tasks_by_name(project, taskname)

        # There are no workflows for the sample so we need to submit one.
        if not tasks:
            logger.debug(
                "No previous %s workflows found for sample %s", workflow_name, sample
            )
            sample_workflow = platform.submit_task(
                taskname, project, workflow, parameters
            )
            if sample_workflow is None:
                logger.error("Failed to submit task for sample %s", sample)
                samples[sample]["workflows"][workflow_name]["state"] = "Failed"
                continue
            if sample_workflow is not None:
                samples[sample]["workflows"][workflow_name]["state"] = "Submitted"
                samples[sample]["workflows"][workflow_name][
                    "workflow"
                ] = sample_workflow
                logger.info("  - Submitted")
                continue

        # Scan through all the workflows for the sample and determine if:
        # 1. If there was a successful run, then consider the sample complete.
        # 2. If all runs of a sample failed, then consider the sample failed.
        # 3. If there is a running sample, then consider sample running
        # 4. If any runs cancelled, ignore it.
        completed_count = 0
        failed_count = 0
        running_count = 0
        cancelled_count = 0
        # Count the state ofall the tasks for the sample to determine the overall state of the sample.
        for task in tasks:
            workflow_state = platform.get_task_state(task)

            if workflow_state == "Complete":
                completed_count += 1
                samples[sample]["workflows"][workflow_name]["workflow"] = task
            elif workflow_state == "Failed":
                failed_count += 1
            elif workflow_state in ["Queued", "Running", "Locked"]:
                running_count += 1
                samples[sample]["workflows"][workflow_name]["workflow"] = task
            elif workflow_state == "Cancelled":
                cancelled_count += 1
            else:
                raise ValueError(f"Unknown workflow state: {workflow_state}")

        # Determine the workflow state of the sample based on the counts above.
        if completed_count > 0:
            samples[sample]["workflows"][workflow_name]["state"] = "Complete"
            logger.info("  - Completed")
            continue

        if running_count > 0:
            samples[sample]["workflows"][workflow_name]["state"] = "Running"
            logger.info("  - Running")
            continue

        if failed_count >= 300:
            samples[sample]["workflows"][workflow_name]["state"] = "Failed"
            logger.info("  - Failed")
            continue

        # If none succeeded, determine if all workflows failed
        if failed_count == len(tasks) - cancelled_count:
            sample_workflow = platform.submit_task(
                taskname, project, workflow, parameters
            )
            if sample_workflow is None:
                logger.error("Failed to submit task for sample %s", sample)
                samples[sample]["workflows"][workflow_name]["state"] = "Failed"
                continue
            samples[sample]["workflows"][workflow_name]["state"] = "Submitted"
            samples[sample]["workflows"][workflow_name]["workflow"] = sample_workflow
            logger.info("  - Submitted")
            continue

        # If the sample was cancelled, and there are no other cases from above, resubmit.
        if cancelled_count == len(tasks):
            sample_workflow = platform.submit_task(
                taskname, project, workflow, parameters
            )
            if sample_workflow is None:
                logger.error("Failed to submit task for sample %s", sample)
                samples[sample]["workflows"][workflow_name]["state"] = "Failed"
                continue
            if sample_workflow is not None:
                samples[sample]["workflows"][workflow_name]["state"] = "Submitted"
                samples[sample]["workflows"][workflow_name][
                    "workflow"
                ] = sample_workflow
                logger.info("  - Submitted")
                continue

    return samples


def wait_for_analysis(platform, samples, workflow_name):
    """
    Wait for a list of tasks to complete
    :param samples: Dictionary of sample name -> container_request_uuid or container_uuid
    :return: Dictionary of sample -> output_uuid
    """
    for counter, sample in enumerate(samples):
        logger.info("[%d/%d] %s", counter + 1, len(samples), sample)

        if workflow_name not in samples[sample]["workflows"]:
            continue

        if "state" not in samples[sample]["workflows"][workflow_name].keys():
            continue

        # If sample is done, skip it.
        if samples[sample]["workflows"][workflow_name]["state"] in (
            "Complete",
            "Cancelled",
            "Failed",
            "NA",
        ):
            continue

        # If not finished, wait for it.
        samples[sample]["workflows"][workflow_name]["state"] = platform.get_task_state(
            samples[sample]["workflows"][workflow_name]["workflow"], refresh=True
        )
        while samples[sample]["workflows"][workflow_name]["state"] not in (
            "Complete",
            "Cancelled",
            "Failed",
        ):
            logger.debug(" - Sleeping for 60 seconds")
            time.sleep(60)
            samples[sample]["workflows"][workflow_name][
                "state"
            ] = platform.get_task_state(
                samples[sample]["workflows"][workflow_name]["workflow"], refresh=True
            )
    return samples


def wait_for_workflow(platform, workflow_name, workflow):
    """
    Wait for a workflow to complete.
    :param workflow: Workflow to wait for
    :return: Workflow state
    """
    logger.info(" - Waiting on %s", workflow_name)
    state = platform.get_task_state(workflow, refresh=True)
    while state not in ("Complete", "Cancelled", "Failed"):
        logger.debug(" - Sleeping for 60 seconds")
        time.sleep(60)
        state = platform.get_task_state(workflow, refresh=True)
    return state


def construct_output_list(platform, samples, merge_workflow):
    """
    Construct a list of files from all the various pipelines to export, where each entry in the list
    is a dictionary with the following keys:
        'source': The source file to copy
        'destination': The destination folder to copy the file to
    """
    # This is also not the right place to hard-code this.
    # This should be in a config file or as tags on the CWL worflow itself.
    outputs = []
    for _, sample in enumerate(samples):
        for workflow, workflow_outputs in pipeline_config.WORKFLOW_TO_OUTPUTS.items():
            if "state" not in samples[sample]["workflows"][workflow].keys():
                continue
            if samples[sample]["workflows"][workflow]["state"] == "Complete":
                for (output,target_folder) in workflow_outputs:
                    outputs.append(
                        {
                            "source": platform.get_task_output(
                                samples[sample]["workflows"][workflow]["workflow"],
                                output,
                            ),
                            "destination": target_folder,
                            "sample": sample,
                            "workflow_output": workflow + "/" + output,
                        }
                    )
    
    merged_output = platform.get_task_output(merge_workflow, "merged_table")
    outputs.append({
        "source": merged_output,
        "destination": "Merged",
        "sample": "Merged",
        "workflow_output": "merge_count_table/merged_table",
    })
    return outputs


def do_work(args, platform, project, launcher_task):
    """Do heavy lifting to run the workflow"""
    # Read samplesheet
    samples = get_samples(args.sample_sheet)

    # Copy reference workflow(s)
    workflows = copy_workflows(
        pipeline_config.PIPELINE_PLATFORM_CONFIG[args.platform], platform, project
    )

    if launcher_task:
        args.config_id = platform.get_task_input(launcher_task, "config")
        args.reffile_id = platform.get_task_input(launcher_task, "reffile")
    else:
        # upload input files to project
        args.config_id = platform.upload_file_to_project(
            args.config,
            project,
            "/launcher_inputs/",
            os.path.basename(args.config),
            overwrite=True,
        )
        args.reffile_id = platform.upload_file_to_project(
            args.reffile,
            project,
            "/launcher_inputs/",
            os.path.basename(args.reffile),
            overwrite=True,
        )

    # Run per-sample  workflow
    parameters = get_default_parameters(args)
    samples = run_persample_workflow(
        samples,
        workflows["per_sample_workflow"],
        parameters,
        "per_sample",
        platform,
        project,
        args.no_reuse
    )
    samples = wait_for_analysis(platform, samples, "per_sample")

    # Run merge count table workflow
    merge_workflow = run_merge_count_table_workflow(
        samples,
        workflows["merge_workflow"],
        platform,
        project,
    )

    if not merge_workflow:
        logger.error("Failed to run merge_metrics_workflow")
    else:
        # Wait for merge_metrics_workflow to complete
        wait_for_workflow(platform, "merge_metrics", merge_workflow)
        rename_output(platform, project, merge_workflow)

    # Stage output files
    outputs = construct_output_list(platform, samples, merge_workflow)

    # Generate list of outputs
    generate_output_file_list(outputs)


def parse_arguments(argv):
    """Parse command line arguments"""
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-l",
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Prints warnings to console by default",
    )
    parser.add_argument(
        "-v", "--version", action="version", version=pipeline_config.__version__
    )
    parser.add_argument("--platform", default=None, choices=SUPPORTED_PLATFORMS.keys())

    # Allow user to provide platform project to run in when running launcher locally.
    # These aren't needed when launcher is run on a platform.
    project = parser.add_mutually_exclusive_group(required=False)
    project.add_argument(
        "--project_name", help="Project name where this workflow is executed"
    )
    project.add_argument(
        "--project_id", help="Project ID/UUID where this workflow is executed"
    )

    parser.add_argument(
        "-s", "--sample_sheet", required=True, help="Samplesheet for pipeline"
    )
    parser.add_argument(
        "--no_reuse", help="Do not reuse results from previous sussessful runs",
        action='store_true', default=False
    )

    parser.add_argument("--reffile", help="Reference sequence file", required=True)

    parser.add_argument("--primer_1", help="Left Primer sequence (5'->3')", required=False)
    parser.add_argument("--primer_2", help="Right Primer sequence (5'->3')", required=False)
    parser.add_argument("--min_qual", type=float, help="Minimal average base quality",
        required=False, default=30)

    parser.add_argument("--config", help="Transgene config file", required=True)
    parser.add_argument("--truncated_cutoff", help="Truncated HDR/ITR Threshold", type=float,
        required=False, default=0.99)
    parser.add_argument("--five_prime_HA_seq", help="5' Homology arm sequence", required=True)
    parser.add_argument("--three_prime_HA_seq", help="3' Homology arm sequence", required=True)
    parser.add_argument("--HA_match_ratio", type=float,
        help="Homology arm sequence match ratio", required=False, default=0.98
    )

    parser.add_argument("--left_itr_seq", type=str,
        help="Left side ITR sequence in 5'->3' direction", required=False)
    parser.add_argument("--right_itr_seq", type=str,
        help="Right side ITR sequence in 5'->3' direction", required=False)
    parser.add_argument("--seed_size", type=int, help="Seed size for blastn", required=False)
    parser.add_argument("--perc_identity", type=int,
        help="Min percent identity to match a ITR sequence for blastn", required=False)

    return parser.parse_args(argv)

def main(argv):
    """main entry point"""
    # Parse arguments
    args = parse_arguments(argv)

    # Start logging
    logger.info("%s v%s", pipeline_config.__appname__, pipeline_config.__version__)
    logger.info(args)

    # Construct and connect to platform
    on_platform = False
    if not args.platform:
        # See if we can figure out what platform we are running on.
        logger.info("No platform provided...detecting platform.")
        args.platform = PlatformFactory().detect_platform()
        logger.info("Detected platform: %s", args.platform)
        on_platform = True

    # Connect to the platform
    platform = PlatformFactory().get_platform(args.platform)
    platform.set_logger(logger)
    platform_config = pipeline_config.PIPELINE_PLATFORM_CONFIG[args.platform]
    platform.connect(**platform_config)

    # Get the project uuid either by its provided name, or from this launcher's project.
    if args.project_id:
        project = platform.get_project_by_id(args.project_id)
    elif args.project_name:
        project = platform.get_project_by_name(args.project_name)
    else:
        project = platform.get_project()

    if not project:
        logger.error("Could not determine project to run in")
        sys.exit(1)

    # If we are running on the platform, get the launcher task
    launcher_task = platform.get_current_task() if on_platform else None

    # do work
    do_work(args, platform, project, launcher_task)

    # done
    logger.info("Done.")


if __name__ == "__main__":
    main(sys.argv[1:])
