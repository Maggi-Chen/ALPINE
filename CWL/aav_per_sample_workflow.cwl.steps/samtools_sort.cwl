cwlVersion: v1.2
class: CommandLineTool
label: Samtools sort index
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: ShellCommandRequirement
- class: ResourceRequirement
  coresMin: 8
  ramMin: 32000
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: aavlr_samtools:1.19.2

inputs:
  samfile:
    type: File
    default: aligned.sam
    inputBinding:
      position: 1
  sample_name:
    type: string
  thread:
    type: int?
    default: 8
    inputBinding:
      position: 2
      prefix: -@

outputs:
  - id: sorted_bam
    type: File
    secondaryFiles:
    - pattern: .bai
      required: false
    outputBinding:
      glob: "$(inputs.sample_name + \".bam\")\n"


baseCommand: []
arguments:
- position: -1
  valueFrom: "samtools sort "
  shellQuote: false
- position: 3
  valueFrom: |-
    ${ 	
       return " -o " + inputs.sample_name + ".bam"
    }
  shellQuote: false
- position: 10
  valueFrom: |-
    ${  
       return " && samtools index " + inputs.sample_name + ".bam"
    }
  shellQuote: false

