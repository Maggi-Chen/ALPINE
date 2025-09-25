cwlVersion: v1.2
class: CommandLineTool
label: Merge Classification tables
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: ShellCommandRequirement
- class: ResourceRequirement
  coresMin: 1
  ramMin: 8000
- class: DockerRequirement
  dockerPull: aavlr_count:1.2
- class: InlineJavascriptRequirement

baseCommand: ["python","/opt/count_classification_table.py"]

inputs:
  - id: classification
    type: File
    label: Classification tables
    doc: Input read classification table for one sample to count reads in each category.
    inputBinding:
      position: 1
      prefix: -i
  - id: output_name
    type: string
    label: Output Filename Prefix
    doc: Output Filename prefix.
    inputBinding:
      position: 5
      prefix: -o

outputs:
  - id: count_table
    type: File
    outputBinding:
      glob: "$(inputs.output_name + \"_read_classification.txt\")\n"
  - id: pie_chart
    type: File
    outputBinding:
      glob: "$(inputs.output_name + \"_read_classification.pdf\")\n"
