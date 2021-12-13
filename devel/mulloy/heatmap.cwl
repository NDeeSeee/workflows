cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement
  - class: ScatterFeatureRequirement


inputs:

  scores_files:
    type: File[]

  scores_labels:
    type: string[]

  regions_files:
    type: File[]

  regions_labels:
    type: string[]

  before_region_start_length:
    type: int
  
  after_region_start_length:
    type: int

  bin_size:
    type: int

  plot_title:
    type: string

  recentering:
    type:
    - "null"
    - type: enum
      symbols: ["Gene TSS", "Peak Center"]
    default: "Gene TSS"

  plot_type:
    type:
    - type: enum
      name: "plot_type"
      symbols: ["lines", "fill", "se", "std"]
  
  sort_regions:
    type:
    - type: enum
      name: "sort_regions"
      symbols: ["descend", "ascend", "no", "keep"]

  what_to_show:
    type:
    - type: enum
      name: "what_to_show"
      symbols:
      - plot, heatmap and colorbar
      - plot and heatmap
      - heatmap only
      - heatmap and colorbar

  per_group:
    type: boolean

  threads:
    type: int


outputs:

  scores_matrix:
    type: File
    outputSource: compute_score_matrix/scores_matrix

  compute_score_matrix_stdout_log:
    type: File
    outputSource: compute_score_matrix/stdout_log

  compute_score_matrix_stderr_log:
    type: File
    outputSource: compute_score_matrix/stderr_log

  make_heatmap_stdout_log:
    type: File
    outputSource: make_heatmap/stdout_log

  make_heatmap_stderr_log:
    type: File
    outputSource: make_heatmap/stderr_log

  heatmap_file:
    type: File
    outputSource: make_heatmap/heatmap_file


steps:

  recenter_regions:
    run: ../../tools/custom-bash.cwl
    in:
      input_file: regions_files
      param: recentering
      script:
        default: |
          if [ "$1" == "Gene TSS" ]
          then
            # BED for gene list
            # chrom  start  end  name  [score] strand
            echo "Recenter by the gene TSS"
            cat "$0" | tr -d "\r" | tr "," "\t" | awk NF | sort -u -k1,1 -k2,2n -k3,3n | awk '{tss=$2; if ($6=="-") tss=$3; print $1"\t"tss"\t"tss+1"\t"$4"\t"$5"\t"$6}' > `basename $0`
          else
            # BED for peaks
            # chrom  start  end
            echo "Recenter by the peak center"
            cat "$0" | tr -d "\r" | tr "," "\t" | awk NF | sort -u -k1,1 -k2,2n -k3,3n | awk '{center=$2+int(($3-$2)/2); print $1"\t"center"\t"center+1}' > `basename $0`
          fi
    scatter: input_file
    out: [output_file]

  compute_score_matrix:
    run: ../../tools/deeptools-computematrix-referencepoint.cwl
    in:
      score_files: scores_files
      regions_files: recenter_regions/output_file
      reference_point: 
        default: "TSS"  # doesn't matter what we set here because we centered regions ourlselves
      before_region_start_length: before_region_start_length
      after_region_start_length: after_region_start_length
      bin_size: bin_size
      sort_regions: sort_regions
      samples_label: scores_labels
      output_filename:
        default: "score_matrix.gz"
      missing_data_as_zero:
        default: true
      threads: threads
    out:
    - scores_matrix
    - stdout_log
    - stderr_log

  make_heatmap:
    run: ../../tools/deeptools-plotheatmap.cwl
    in:
      plot_title: plot_title
      scores_matrix: compute_score_matrix/scores_matrix
      output_filename:
        default: "score_matrix.pdf"
      plot_type: plot_type
      sort_regions: sort_regions
      average_type_summary_plot:
        default: "mean"
      what_to_show: what_to_show
      ref_point_label: recentering
      regions_label: regions_labels
      samples_label: scores_labels
      y_axisLabel:
        default: "Signal mean"
      per_group: per_group
      legend_location:
        default: "upper-left"
    out:
    - heatmap_file
    - stdout_log
    - stderr_log
