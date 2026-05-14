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
    inputBinding:
      position: 3
      prefix: -x
      valueFrom: |
        ${
          if (inputs.preset !== null && inputs.preset !== undefined) {
            return inputs.preset;
          } else {
            var data_type = inputs.data_type;
            if (data_type === "pacbio-hifi") {
              return "map-hifi";
            } else if (data_type === "nanopore") {
              return "map-ont";
            } else {
              return "map-hifi";
            }
          }
        }
  thread:
    type: int?
    default: 8
    inputBinding:
      position: 4
      prefix: -t
  data_type:
    type:
      - 'null'
      - name: data_choices
        type: enum
        symbols:
          - pacbio-hifi
          - nanopore
    default: "pacbio-hifi"


outputs:
  - id: aligned_sam
    type: File
    outputBinding: 
      glob: aligned.sam


