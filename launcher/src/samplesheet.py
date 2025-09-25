""" Samplesheet functions """
import logging

logger = logging.getLogger(__name__)

def generate_output_file_list(outputs, output_file="output_files.txt"):
    """
    Generate list of output files
    # format is:
    # sampleid     workflow/output_field     output_file_id
    """
    with open(output_file, mode="w", encoding="utf-8") as file:
        for output in outputs:
            file.write(
                "\t".join(
                    [output["sample"], output["workflow_output"], output["source"]]
                )
                + "\n"
            )


def get_samples(sample_sheet_file):
    """
    Read the sample sheet and return a samples dictionary
    Input format: 2 fields:
    sample_id, fastq_file

    Return a dictionary of all samples, where key is sample_id, e.g.
    sample -> {
            'fastq' = ...
        }
    """
    # read samplesheet and populate sample to subject map
    # create subject to sample map from samplesheet rows
    samples = {}
    with open(sample_sheet_file, "r", encoding="UTF-8") as file:
        for line in file.readlines():
            line = line.strip()
            if not line:
                continue
            (
                sample,
                fastq_file,
            ) = line.split("\t")

            samples[sample] = {}
            samples[sample]['fastq'] = fastq_file

    return samples
