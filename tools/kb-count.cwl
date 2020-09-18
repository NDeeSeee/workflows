cwlVersion: v1.0
class: CommandLineTool


hints:
  - class: DockerRequirement
    dockerPull: biowardrobe2/kb-python:v0.0.1


inputs:

  fastq_file_1:
    type: File
    inputBinding:
      position: 100
    doc: "Fastq file 1"

  fastq_file_2:
    type: File
    inputBinding:
      position: 101
    doc: "Fastq file 2"
  
  fastq_file_3:
    type: File?
    inputBinding:
      position: 102
    doc: "Fastq file 3"
    
  kallisto_index_file:
    type: File
    inputBinding:
      position: 5
      prefix: "-i"
    doc: "Kallisto index file"

  transcript_to_gene_mapping_file:
    type: File
    inputBinding:
      position: 6
      prefix: "-g"
    doc: "Transcript-to-gene mapping file"

  sc_technology:
    type:
    - type: enum
      name: "sc_technology"
      symbols:
      - 10XV1       # 3 input files
      - 10XV2       # 2 input files 
      - 10XV3       # 2 input files 
      - CELSEQ      # 2 input files
      - CELSEQ2     # 2 input files
      - DROPSEQ     # 2 input files
      - INDROPSV1   # 2 input files
      - INDROPSV2   # 2 input files
      - INDROPSV3   # 3 input files
      - SCRUBSEQ    # 2 input files
      - SURECELL    # 2 input files
    inputBinding:
      position: 7
      prefix: "-x"
    doc: "Single-cell technology used"

  loom:
    type: boolean?
    inputBinding:
      position: 8
      prefix: "--loom"
    doc: "Generate loom file from count matrix"

  h5ad:
    type: boolean?
    inputBinding:
      position: 9
      prefix: "--h5ad"
    doc: "Generate h5ad file from count matrix"


outputs:

  counts_unfiltered_folder:
    type: Directory
    outputBinding:
      glob: "counts_unfiltered"
    doc: |
      Count matrix files generated by bustools count
      command

  whitelist_file:
    type: File
    outputBinding:
      glob: "*_whitelist.txt"
    doc: |
      File of whitelisted barcodes. Corresponds to the used
      single-cell technology

  bustools_inspect_report:
    type: File
    outputBinding:
      glob: "inspect.json"
    doc: |
      Report summarizing the contents of a sorted BUS file.
      Result of running bustools sort with not_sorted_bus_file,
      then bustools inspect with sorted BUS file

  kallisto_bus_report:
    type: File
    outputBinding:
      glob: "run_info.json"
    doc: |
      Report generated by kallisto bus run

  ec_mapping_file:
    type: File
    outputBinding:
      glob: "matrix.ec"
    doc: |
      File for mapping equivalence classes to transcripts.
      Direct output of kallisto bus command

  transcripts_file:
    type: File
    outputBinding:
      glob: "transcripts.txt"
    doc: |
      File to store transcript names.
      Direct output of kallisto bus command

  not_sorted_bus_file:
    type: File
    outputBinding:
      glob: "output.bus"
    doc: |
      Not sorted BUS file.
      Direct output of kallisto bus command

  corrected_sorted_bus_file:
    type: File
    outputBinding:
      glob: "output.unfiltered.bus"
    doc: |
      Sorted BUS file with corrected barcodes.
      Result of running bustools sort with
      not_sorted_bus_file, then bustools correct
      using whitelist_file and bustools sort with
      corrected BUS file


baseCommand: ["kb", "count", "--verbose"]


doc: |
  Uses kallisto to pseudoalign reads and bustools to quantify the data.

  1. Generates BUS file from input fastq files
  2. Sorts generated BUS file
  3. Inspects sorted BUS file
  4. Corrects barcodes in sorted BUS file
  5. Sorts corrected BUS file
  6. Generates count matrix from sorted barcode corrected BUS file

  Notes:
  --verbose was hardcoded
  --lamanno and --nucleus arguments were skipped, so we don't need -c1, -c2
  --keep-tmp, --overwrite doesn't make sense when running from container