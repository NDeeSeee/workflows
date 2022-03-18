cwlVersion: v1.1
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing: |
    ${
      return [
        {
          "entry": inputs.annotation_gtf_file,
          "entryname": inputs.annotation_gtf_file.basename,
          "writable": true
        }
      ]
    }


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/spladder:v0.0.1


inputs:

  alignment_files:
    type:
    - File
    - File[]
    secondaryFiles:
    - .bai
    inputBinding:
      prefix: "--bams"
    doc: |
      Sorted and indexed alignment file(s) in BAM format

  annotation_gtf_file:
    type: File
    inputBinding:
      prefix: "--annotation"
    doc: |
      Annotation file in GTF or GFF3 format

  merge_strategy:
    type:
    - "null"
    - type: enum
      symbols:
      - "single"
      - "merge_bams"
      - "merge_graphs"
      - "merge_all"
    inputBinding:
      prefix: "--merge-strat"
    doc: |
      Merge strategy if multiple BAM files are provided.
      single - generates a separate augmented splicing graph per given input alignment files.
      merge_bams - treats all input BAM files as technical replicates and directly forms a splicing graph using all reads.
      merge_graphs - generates a single splicing graph per input BAM file and then merges all graphs into a joint single graph.
      merge_all - a combination of merge_bams and merge_graphs modes.
      Default: merge_graphs

  ignore_mismatches:
    type: boolean?
    inputBinding:
      prefix: "--ignore-mismatches"
    doc: |
      Do not filter by edit operations - do not require NM tag in BAM.
      Default: SplAdder expects the NM tag in all BAM files

  threads:
    type: int?
    inputBinding:
      prefix: "--parallel"
    doc: |
      Number of cores/cpus to use.
      Default: 1


outputs:

  # results_folder:
  #   type: Directory
  #   outputBinding:
  #     glob: "./"

  alt_3prime_C3_counts_hdf5:
    type: File?
    outputBinding:
      glob: "./*_alt_3prime_C3.counts.hdf5"

  alt_5prime_C3_counts_hdf5:
    type: File?
    outputBinding:
      glob: "./*_alt_5prime_C3.counts.hdf5"

  exon_skip_C3_counts_hdf5:
    type: File?
    outputBinding:
      glob: "./*_[!mult]_exon_skip_C3.counts.hdf5"

  intron_retention_C3_counts_hdf5:
    type: File?
    outputBinding:
      glob: "./*_intron_retention_C3.counts.hdf5"

  mult_exon_skip_C3_counts_hdf5:
    type: File?
    outputBinding:
      glob: "./*_mult_exon_skip_C3.counts.hdf5"

  mutex_exons_C3_counts_hdf5:
    type: File?
    outputBinding:
      glob: "./*_mutex_exons_C3.counts.hdf5"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["spladder", "build", "--outdir", "./"]


stdout: spladder_build_stdout.log
stderr: spladder_build_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "SplAdder Build"
s:name: "SplAdder Build"
s:alternateName: "Constructs splicing graph and extracts alternative splicing events"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/spladder-build.cwl
s:codeRepository: https://github.com/Barski-lab/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Cincinnati Children's Hospital Medical Center"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: "45229"
    s:streetAddress: "3333 Burnet Ave"
    s:telephone: "+1(513)636-4200"
  s:logo: "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
  s:department:
  - class: s:Organization
    s:legalName: "Allergy and Immunology"
    s:department:
    - class: s:Organization
      s:legalName: "Barski Research Lab"
      s:member:
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898


doc: |
  SplAdder Build
  ===================================================================
  Constructs splicing graph and extracts alternative splicing events.


s:about: |
  usage: spladder build [-h] -b FILE1,FILE2,... -o DIR -a FILE [--parallel <INT>] [-v] [-d] [-n INT] [--primary-only] [--no-primary-only] [--var-aware] [--no-var-aware]
                        [--set-mm-tag MM_TAG] [--labels STRING] [--filter-overlap-genes] [--filter-overlap-exons] [--filter-overlap-transcripts] [--filter-consensus STRING]
                        [--ignore-mismatches] [--reference FILE] [-l FILE] [--output-txt] [--output-txt-conf] [--no-output-txt-conf] [--output-gff3] [--output-gff3-conf]
                        [--no-output-gff3-conf] [--output-struc] [--output-struc-conf] [--output-bed] [--output-conf-bed] [--output-conf-tcga] [--output-conf-icgc] [--sparse-bam]
                        [--compress-text] [--no-compress-text] [--tmp-dir DIR] [-c INT] [-I INT] [-M <STRAT>]
                        [--chunked-merge LEVEL MAX_LEVEL START END [LEVEL MAX_LEVEL START END ...]] [--chunksize INT] [--insert-ir] [--no-insert-ir] [--insert-es] [--no-insert-es]
                        [--insert-ni] [--no-insert-ni] [--remove-se] [--no-remove-se] [--validate-sg] [--no-validate-sg] [--re-infer-sg] [--no-re-infer-sg]
                        [--validate-sg-count INT] [--event-types STRING] [--extract-ase] [--no-extract-ase] [--ase-edge-limit INT] [--curate-alt-prime] [--no-curate-alt-prime]
                        [--quantify-graph] [--no-quantify-graph] [--use-anno-support] [--no-use-anno-support] [--psi-min-reads INT] [--qmode STRING]

  optional arguments:
    -h, --help            show this help message and exit

  MANDATORY:
    -b FILE1,FILE2,..., --bams FILE1,FILE2,...
                          alignment files in BAM format (comma separated list)
    -o DIR, --outdir DIR  output directory
    -a FILE, --annotation FILE
                          file name for annotation in GTF/GFF3 or format

  GENERAL:
    --parallel <INT>      use INT processors [1]
    -v, --verbose         use verbose output mode [off]
    -d, --debug           use debug output mode [off]

  INPUT OPTIONS:
    -n INT, --readlen INT
                          read length (used for automatic confidence levele settings) [36]
    --primary-only        only use primary alignments [on]
    --no-primary-only
    --var-aware           alignment files are variation aware (presence of XM and XG tags) [off]
    --no-var-aware
    --set-mm-tag MM_TAG   sets the sequence of the mismatch tag used in alignments [NM]
    --labels STRING       use labels instead of bam file names (comma separated list) [-]
    --filter-overlap-genes
                          remove genes from annotation that overlap each other [off]
    --filter-overlap-exons
                          remove exons from annotation that overlap each other [off]
    --filter-overlap-transcripts
                          remove transcripts from annotation that overlap each other [off]
    --filter-consensus STRING
                          require new junctions to have consensus (needs ref genome) [off]; choices: strict (GT/AG), lenient (G[TC]/AG)
    --ignore-mismatches   does not filter by edit operations - does not require NM in BAM [off]
    --reference FILE      reference genome (only needed for CRAM file de-compression or consensus filtering)

  OUTPUT OPTIONS:
    -l FILE, --logfile FILE
                          log file name [stdout]
    --output-txt          outputs all events in txt format (can be big) [off]
    --output-txt-conf     outputs confirmed events in txt format [on]
    --no-output-txt-conf
    --output-gff3         outputs all events in GFF3 format [off]
    --output-gff3-conf    outputs confirmed events in GFF3 format [on]
    --no-output-gff3-conf
    --output-struc        outputs all events in structured splicing syntax similar to astalavista [off]
    --output-struc-conf   outputs confirmed events in structured splicing syntax similar to astalavista [off]
    --output-bed          output all events in BED format [off]
    --output-conf-bed     output confirmed events in BED format [off]
    --output-conf-tcga    output confirmed events in format used for TCGA [off]
    --output-conf-icgc    output confirmed events in format used for ICGC [off]
    --sparse-bam          store BAM content as sparse representation for later use [off]
    --compress-text       compress text output [on]
    --no-compress-text
    --tmp-dir DIR         directory to store temporary data [<outdir>/tmp]

  GRAPH GENERATION:
    -c INT, --confidence INT
                          confidence level (0 lowest to 3 highest) [3]
    -I INT, --iterations INT
                          number of iterations to insert new introns into the graph [5]
    -M <STRAT>, --merge-strat <STRAT>
                          merge strategy, where <STRAT> is one of: single, merge_bams, merge_graphs, merge_all [merge_graphs]
    --chunked-merge LEVEL MAX_LEVEL START END [LEVEL MAX_LEVEL START END ...]
                          provide info for external merge with START being 0-based and END non-inclusive
    --chunksize INT       chunksize for chunked merge [10]
    --insert-ir           insert intron retentions [on]
    --no-insert-ir
    --insert-es           insert cassette exons [on]
    --no-insert-es
    --insert-ni           insert new intron edges [on]
    --no-insert-ni
    --remove-se           remove short exons [off]
    --no-remove-se
    --validate-sg         validate splice graph [off]
    --no-validate-sg
    --re-infer-sg         re-infer splice graph [off] (advanced)
    --no-re-infer-sg
    --validate-sg-count INT
                          number of samples supporting an edge for it to be kept [min(10, #samples)]

  AS EVENT EXTRACTION:
    --event-types STRING  list of alternative splicing events to extract [exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip,mutex_exons]
    --extract-ase         extract alternative splicing events [on]
    --no-extract-ase
    --ase-edge-limit INT  max number of edges in the graph to still extract events for a gene [500]
    --curate-alt-prime    curate alt prime events [on]
    --no-curate-alt-prime
    --quantify-graph      quantify graph [on]
    --no-quantify-graph
    --use-anno-support    use annotation for validating event introns [off]
    --no-use-anno-support
    --psi-min-reads INT   minimum number of spliced reads covering either isoform to compute PSI [10]
    --qmode STRING        quantification mode: single, collect, all [all]