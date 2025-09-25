cwlVersion: v1.2
class: CommandLineTool
label: Merge Classification tables v1.4
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: ShellCommandRequirement
- class: ResourceRequirement
  coresMin: 1
  ramMin: 8000
- class: DockerRequirement
  dockerPull: aavlr_merge:1.2
- class: InlineJavascriptRequirement

baseCommand: ["python","/opt/merge_count_table.py"]

inputs:
  - id: count_tables
    type: File[]
    label: Count tables
    doc: Input read count tables for individual samples to merge into one table.
    inputBinding:
      position: 1
      prefix: -i
  - id: output_name
    type: string?
    label: Output Filename Prefix
    doc: Output Filename prefix.
    default: Merged
    inputBinding:
      position: 5
      prefix: -o

outputs:
  - id: merged_table
    type: File
    outputBinding:
      glob: "$(inputs.output_name + \"_read_classification.txt\")\n"
  - id: pie_chart
    type: File
    outputBinding:
      glob: "$(inputs.output_name + \"_read_classification.pdf\")\n"
