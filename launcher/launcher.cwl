cwlVersion: v1.2
class: CommandLineTool
label: Transgene Integration Long-Read Launcher
$namespaces:
  arv: "http://arvados.org/cwl#"

hints:
  DockerRequirement:
    dockerPull: transgene-launcher:__VERSION__

requirements:
  MultipleInputFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
  WorkReuse:
    enableReuse: false
  arv:APIRequirement: {}

baseCommand: ["/opt/launcher/launcher.py"]

inputs:
  sample_sheet:
    type: File
    label: Sample Sheet
    doc: A tab-delimited list of samples containing two columns - Sample_ID and FASTQ_Path.
    inputBinding:
      prefix: -s
      shellQuote: false
  reffile:
    type: File
    label: Reference Sequence File
    doc: Reference sequence file should contain possible sequences for WT (Unmodified), HDR, ITR integration. 
    inputBinding:
      prefix: --reffile
      shellQuote: false
  config:
    type: File
    label: Transgene Config File
    doc: Config file listing detailed information for all reference sequences provided in the Reference file. Should include 9 columns, including AAV_Vector, Ref_Name, Ref_Type, ITR1_start, ITR1_end, gene_start, gene_end, ITR2_start, ITR2_end. A header line is required for the config file.
    inputBinding:
      prefix: --config
      shellQuote: false
  primer_1:
    type: string?
    label: Forward Primer Sequence
    doc: Forward primer sequence in 5'->3' direction.
    inputBinding:
      prefix: --primer_1
      shellQuote: false
  primer_2:
    type: string?
    label: Reverse Primer Sequence
    doc: Reverse primer sequence in 5'->3' direction.
    inputBinding:
      prefix: --primer_2
      shellQuote: false
  min_qual:
    type: float?
    label: Minimal average base quality
    doc: Sequencing reads with average Phred base quality below this cutoff will be filtered out. Default is 30.
    inputBinding:
      prefix: --min_qual
      shellQuote: false
  five_prime_HA_seq:
    type: string
    label: 5' Homology Arm Sequence
    doc: 5' Homology Arm sequence in 5'->3' direction.
    inputBinding:
      prefix: --five_prime_HA_seq
      shellQuote: false
  three_prime_HA_seq:
    type: string
    label: 3' Homology Arm Sequence
    doc: 3' Homology Arm sequence in 5'->3' direction.
    inputBinding:
      prefix: --three_prime_HA_seq
      shellQuote: false
  truncated_cutoff:
    type: float?
    label: Truncated HDR/ITR Threshold
    doc: Truncated HDR/ITR threshold. Should be 0-1, default is 0.95.
    inputBinding:
      prefix: --truncated_cutoff
      shellQuote: false
  HA_Seq_match_ratio:
    type: float?
    label: Homology Arm Sequence Match Ratio
    doc: Homology Arm Sequence Match Ratio. Default is 0.96.
    inputBinding:
      prefix: --HA_match_ratio
      shellQuote: false
  left_itr_seq:
    type: string?
    label: Left Side ITR Sequence
    doc: Left ITR sequence in 5'->3' direction.
    inputBinding:
      prefix: --left_itr_seq
      shellQuote: false
  right_itr_seq:
    type: string?
    label: Right Side ITR Sequence
    doc: Right ITR sequence in 5'->3' direction.
    inputBinding:
      prefix: --right_itr_seq
      shellQuote: false
  disable_reuse:
    type: boolean?
    label: Disable Task Reuse
    doc: Disable task reuse and force launcher to re-process all samples.
    default: False
    inputBinding:
      prefix: --no_reuse
      shellQuote: false

outputs:
  output_files:
    type: File
    outputBinding:
      glob: output_files.txt
