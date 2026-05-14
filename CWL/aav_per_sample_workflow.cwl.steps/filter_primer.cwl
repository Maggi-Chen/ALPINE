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
  dockerPull: aavlr_filter:1.6
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
  - id: primer_check_length
    type: int?
    inputBinding:
      position: 5
      prefix: --primer_check_length
  - id: data_type
    type:
      - 'null'
      - name: data_choices
        type: enum
        symbols:
          - pacbio-hifi
          - nanopore
    default: "pacbio-hifi"

arguments:
- position: 6
  valueFrom: |-
    ${
        if (inputs.min_qual !== null && inputs.min_qual !== undefined) {
            return "  --min_qual " + inputs.min_qual;
        } else {
            if (inputs.data_type == "nanopore") {
                return "  --min_qual 10 ";
            } else {
                return " ";
            }
        }
    }
  shellQuote: false


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
