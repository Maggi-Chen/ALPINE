class: CommandLineTool
cwlVersion: v1.2
$namespaces:
  sbg: 'https://sevenbridges.com'
baseCommand: []
inputs:
  - id: in_array
    type: File[]?
    label: Input array
    doc: Array of input files.
  - id: output_name
    type: string?
    label: Output name
    doc: Name for the output directory
outputs:
  - id: outputs
    doc: Output directory containing all files from input array
    label: Output directory
    type: Directory
    outputBinding:
      glob: $(inputs.output_name)
      loadListing: deep_listing
label: Group Outputs
arguments:
  - prefix: ''
    shellQuote: false
    position: 2
    valueFrom: |-
      ${
          var out_folder = inputs.output_name
          var cmd = "&& while read p; do cp $p " + out_folder +"; done < paths.txt"
          return ["mkdir", out_folder, cmd].join(" ")
      }
requirements:
  - class: ShellCommandRequirement
  - class: LoadListingRequirement
  - class: DockerRequirement
    dockerPull: ubuntu:16.04
  - class: InitialWorkDirRequirement
    listing:
      - entryname: paths.txt
        entry: |-
          ${
              return inputs.in_array.filter(function(a) { return a !== null && a !== undefined; }).map(function(a) {return a.path}).join("\n") + "\n"
          }
        writable: false
  - class: InlineJavascriptRequirement
