cwlVersion: v1.2
class: CommandLineTool
label: Classify reads
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: ShellCommandRequirement
- class: ResourceRequirement
  coresMin: 8
  ramMin: 32000
- class: DockerRequirement
  dockerPull:  aavlr_classify:1.4.1
- class: InlineJavascriptRequirement

baseCommand: ["python","/opt/classify_read.py"]

inputs:
  - id: bamfile
    type: File
    secondaryFiles:
    - pattern: .bai
      required: false
    inputBinding:
      position: 1
      prefix: -i
  - id: reffile
    type: File
    inputBinding:
      position: 2
      prefix: -r
  - id: config
    type: File
    inputBinding:
      position: 3
      prefix: --config
  - id: sample_name
    type: string
    inputBinding:
      position: 4
      prefix: -p
  - id: truncated_cutoff
    type: float?
    inputBinding:
      position: 6
      prefix: --full_hdr_cutoff
  - id: filtered_fastq
    type: File
    inputBinding:
      position: 5
      prefix: -fq
  - id: five_prime_HA_seq
    type: string
    inputBinding:
      position: 7
      prefix: --five_prime_HA_arm_seq
  - id: three_prime_HA_seq
    type: string
    inputBinding:
      position: 8
      prefix: --three_prime_HA_arm_seq
  - id: HA_Seq_match_ratio
    type: float?
    inputBinding:
      position: 9
      prefix: --ha_match_ratio
  - id: left_itr_seq
    type: string?
    inputBinding:
      position: 10
      prefix: --left_itr_seq
  - id: right_itr_seq
    type: string?
    inputBinding:
      position: 11
      prefix: --right_itr_seq
  - id: seed_size
    type: int?
    inputBinding:
      position: 12
      prefix: --seed_size
  - id: perc_identity
    type: int?
    inputBinding:
      position: 13
      prefix: --perc_identity
  
outputs:
  - id: read_classification
    type: File
    outputBinding:
      glob: "$(\"readname_\"+inputs.sample_name+\".txt\")\n"
