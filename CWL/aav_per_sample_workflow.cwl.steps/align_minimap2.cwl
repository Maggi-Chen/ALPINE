cwlVersion: v1.2
class: CommandLineTool
label: Minimap2
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: ShellCommandRequirement
- class: ResourceRequirement
  coresMin: 8
  ramMin: 32000
- class: DockerRequirement
  dockerPull: aavlr_minimap2:2.24

baseCommand: ["minimap2","-a","--eqx","-o","aligned.sam"]

inputs:
  reffile:
    type: File
    inputBinding:
      position: 0
  fastqfile:
    type: File
    inputBinding:
      position: 1
  preset:
    type: string?
    default: "map-hifi"
    inputBinding:
      position: 3
      prefix: -x
  thread:
    type: int?
    default: 8
    inputBinding:
      position: 4
      prefix: -t


outputs:
  - id: aligned_sam
    type: File
    outputBinding: 
      glob: aligned.sam


