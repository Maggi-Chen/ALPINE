cwlVersion: v1.2
class: CommandLineTool
label: Filter fastq
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: ShellCommandRequirement
- class: ResourceRequirement
  coresMin: 8
  ramMin: 32000
- class: DockerRequirement
  dockerPull: aavlr_filter:1.2
- class: InlineJavascriptRequirement

baseCommand: ["python","/opt/filter_fastq_gz.py"]

inputs:
  - id: fastqfile
    type: File
    inputBinding:
      position: 1
      prefix: -i
  - id: primer_f
    type: string?
    inputBinding:
      position: 2
      prefix: -f
  - id: primer_r
    type: string?
    inputBinding:
      position: 3
      prefix: -r
  - id: sample_name
    type: string
    inputBinding:
      position: 4
      prefix: --sample
  - id: min_qual
    type: float?
    inputBinding:
      position: 6
      prefix: --min_qual

outputs:
  - id: filtered_fastq
    type: File
    outputBinding:
      glob: "$(inputs.sample_name + \".fastq.gz\")\n"
  - id: filter_log
    type: File
    outputBinding:
      glob: "$(inputs.sample_name + \".filtering.log\")\n"
  - id: original_histogram_plot
    type: File
    outputBinding:
      glob: "*_original_combined_histograms.html"
  - id: filtered_histogram_plot
    type: File
    outputBinding:
      glob: "*_filtered_combined_histograms.html"
