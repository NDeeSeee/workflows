cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scvelo:v0.0.1


inputs:

  adata_file:
    type:
    - File
    - File[]
    inputBinding:
      prefix: "--h5ad"
    doc: |
      Path to the h5ad file.
      If multiple files provided, the order should correspond to
      the order of files provided in --cells, --umap, --pca, --clusters

  cells_metadata:
    type:
    - File
    - File[]
    inputBinding:
      prefix: "--cells"
    doc: |
      Path to the cells data TSV file.
      If multiple files provided, the order should correspond to
      the order of files provided in --h5ad, --umap, --pca, --clusters

  umap_metadata:
    type:
    - File
    - File[]
    inputBinding:
      prefix: "--umap"
    doc: |
      Path to the custom UMAP data TSV file.
      If multiple files provided, the order should correspond to
      the order of files provided in --cells, --h5ad, --pca, --clusters

  pca_metadata:
    type:
    - File
    - File[]
    inputBinding:
      prefix: "--pca"
    doc: |
      Path to the custom PCA data TSV file.
      If multiple files provided, the order should correspond to
      the order of files provided in --cells, --umap, --h5ad, --clusters

  clusters_metadata:
    type:
    - File
    - File[]
    inputBinding:
      prefix: "--clusters"
    doc: |
      Path to the clusters data TSV file.
      If multiple files provided, the order should correspond to
      the order of files provided in --cells, --umap, --pca, --h5ad

  nneighbors:
    type: int?
    inputBinding:
      prefix: "--nneighbors"
    doc: |
      Number of nearest neighbors to use when computing moments for
      velocity estimation
      Default: 30

  top_n_genes:
    type: int?
    inputBinding:
      prefix: "--top"
    doc: |
      Top N genes to display on plots and return in tables
      Default: 10

  selected_genes:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--genes"
    doc: |
      List of genes of interest.
      Default: None

  threads:
    type: int?
    inputBinding:
      prefix: "--threads"
    doc: |
      Number of jobs to run in parallel
      Default: 1

  dpi:
    type: int?
    inputBinding:
      prefix: "--dpi"
    doc: |
      Resolution for exported plots
      Default: 300

  output_prefix:
    type: string?
    inputBinding:
      prefix: "--output"
    doc: |
      Prefix for generated output files.
      Default: ./scvelo_


outputs:

  compressed_adata_file:
    type: File
    outputBinding:
      glob: "*compressed.h5ad"

  putative_driver_genes:
    type: File
    outputBinding:
      glob: "*putative_driver_genes.tsv"

  driver_genes_heatmap:
    type: File?
    outputBinding:
      glob: "*heatmap_driver_genes.png"

  latent_time_plot:
    type: File?
    outputBinding:
      glob: "*latent_time.png"

  paga_plot:
    type: File?
    outputBinding:
      glob: "*paga.png"

  proportions_chart:
    type: File?
    outputBinding:
      glob: "*proportions_chart.png"

  velocity_confidence_plot:
    type: File?
    outputBinding:
      glob: "*velocity_confidence.png"

  velocity_grid_plot:
    type: File?
    outputBinding:
      glob: "*velocity_grid.png"

  velocity_length_plot:
    type: File?
    outputBinding:
      glob: "*velocity_length.png"

  velocity_metrics_plot:
    type: File?
    outputBinding:
      glob: "*velocity_metrics.png"

  velocity_stream_plot:
    type: File?
    outputBinding:
      glob: "*velocity_stream.png"

  gene_expression_plots:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*expression.png"

  gene_phase_plots:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*phase.png"

  selected_gene_expression_plots:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*expression_selected.png"

  selected_gene_phase_plots:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*phase_selected.png"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["run_scvelo.py"]

stdout: run_scvelo_stdout.log
stderr: run_scvelo_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "scVelo - RNA velocity generalized through dynamical modeling"
s:name: "scVelo - RNA velocity generalized through dynamical modeling"
s:alternateName: "scVelo is a scalable toolkit for RNA velocity analysis in single cells, based on Bergen et al. (Nature Biotech, 2020)"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/scvelo.cwl
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
  scVelo - RNA velocity generalized through dynamical modeling
  ============================================================
  scVelo is a scalable toolkit for RNA velocity analysis in single cells,
  based on Bergen et al. (Nature Biotech, 2020)


s:about: |
  usage: run_scvelo.py [-h] --h5ad [H5AD [H5AD ...]] --cells [CELLS [CELLS ...]] --umap [UMAP [UMAP ...]]
                            --pca [PCA [PCA ...]] --clusters [CLUSTERS [CLUSTERS ...]] [--nneighbors NNEIGHBORS]
                            [--genes [GENES [GENES ...]]] [--dpi DPI] [--threads THREADS] [--top TOP] [--output OUTPUT]

  optional arguments:
    -h, --help            show this help message and exit
    --h5ad [H5AD [H5AD ...]]
                          Path to the h5ad file. If multiple files provided, the order should correspond to the order of
                          files provided in --cells, --umap, --pca, --clusters
    --cells [CELLS [CELLS ...]]
                          Path to the cells data TSV file. If multiple files provided, the order should correspond to the
                          order of files provided in --h5ad, --umap, --pca,
                          --clusters
    --umap [UMAP [UMAP ...]]
                          Path to the custom UMAP data TSV file. If multiple files provided, the order should correspond
                          to the order of files provided in --cells, --h5ad, --pca,
                          --clusters
    --pca [PCA [PCA ...]]
                          Path to the custom PCA data TSV file. If multiple files provided, the order should correspond to
                          the order of files provided in --cells, --umap, --h5ad,
                          --clusters
    --clusters [CLUSTERS [CLUSTERS ...]]
                          Path to the clusters data TSV file. If multiple files provided, the order should correspond to the
                          order of files provided in --cells, --umap, --pca,
                          --h5ad
    --nneighbors NNEIGHBORS
                          Number of nearest neighbors to use when computing moments for velocity estimation
    --genes [GENES [GENES ...]]
                          List of genes of interest
    --dpi DPI             Resolution for exported plots
    --threads THREADS     Number of jobs to run in parallel
    --top TOP             Top N genes to display on plots and return in tables
    --output OUTPUT       Prefix for generated output files
