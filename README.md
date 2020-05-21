A simple tutorial for a complex ComplexHeatmap
================
Kevin Blighe
2020-05-21

# 1, introduction

The data used for the heatmap is taken from a large bulk RNA-seq study
on rheumatoid arthritis synovial biopsies, a project on which I co-led
as Senior Bioinformatician. The work was published as [Molecular
Portraits of Early Rheumatoid Arthritis Identify Clinical and Treatment
Response Phenotypes](https://pubmed.ncbi.nlm.nih.gov/31461658/) (Lewis
et al. (2019)).

*ComplexHeatmap* (Gu, Eils, and Schlesner (2016)) is an R Programming
Language (R Core Team (2020)) package that is currently listed in the
[Bioconductor](https://bioconductor.org/) package repository.

# 2, install and load required packages

``` r
  require(RColorBrewer)
  require(ComplexHeatmap)
  require(circlize)
  require(digest)
  require(cluster)
```

If all load successfully, proceed to **Part 3**. Otherwise, go through
the following code chunks in order to ensure that each package is
installed and loaded properly.

*BiocManager* (Morgan (2019))

``` r
  if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager')
```

*RColorBrewer* (Neuwirth (2014))

``` r
  if (!requireNamespace('RColorBrewer', quietly = TRUE))
    BiocManager::install('RColorBrewer')

  require(RColorBrewer)
```

*ComplexHeatmap* (Gu, Eils, and Schlesner (2016))

``` r
  if (!requireNamespace('ComplexHeatmap', quietly = TRUE))
    BiocManager::install('ComplexHeatmap')

  require(ComplexHeatmap)
```

*circlize* (Gu et al. (2014))

``` r
  if (!requireNamespace('circlize', quietly = TRUE))
    BiocManager::install('circlize')

  require(circlize)
```

*digest* (Antoine Lucas et al. (2020))

``` r
  if (!requireNamespace('digest', quietly = TRUE))
    BiocManager::install('digest')

  require(digest)
```

*cluster* (Maechler et al. (2019))

``` r
  if (!requireNamespace('cluster', quietly = TRUE))
    BiocManager::install('cluster')

  require(cluster)
```

# 3, obtain the data

Sample information and raw data FASTQ files are stored in the public
domain under accessions
[E-MTAB-6141](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6141/)
and [PRJEB23131](https://www.ebi.ac.uk/ena/data/view/PRJEB23131),
respectively.

I have separately stored an expression matrix and sample metadata in my
GitHub repository
([github.com/kevinblighe/](https://github.com/kevinblighe/)). The data
will be downloaded to temporary files that will later be deleted by your
operating system. The files are presented as uncompressed, plain text in
order to ensure compatibility across all platforms outside of a standard
build system.

``` r
  tmpfile <- tempfile()
  download.file('https://github.com/kevinblighe/E-MTAB-6141/raw/master/rdata/mat.tsv',
    tmpfile, method = 'auto')
  mat <- read.table(tmpfile, sep = '\t', row.names = 1,
    header = TRUE, stringsAsFactors = FALSE)

  tmpfile <- tempfile()
  download.file('https://github.com/kevinblighe/E-MTAB-6141/raw/master/rdata/metadata.tsv',
    tmpfile, method = 'auto')
  metadata <- read.table(tmpfile, sep = '\t', row.names = 1,
    header = TRUE, stringsAsFactors = FALSE)

  tmpfile <- tempfile()
  download.file('https://github.com/kevinblighe/E-MTAB-6141/raw/master/rdata/sig_genes.list',
    tmpfile, method = 'auto')
  sig_genes <- read.table(tmpfile, sep = '\t',
    header = FALSE, stringsAsFactors = FALSE)[,1]
```

Check the md5 checksums to ensure data integrity / security. The
checksums should be:

  - ‘mat’ object: 1fbbe9568738577a2f3e3dc42e6c75cf
  - ‘metadata’ object: 542a40ccf8b14c51ffa45361c5d3aed9
  - ‘sig\_genes’ object: fdc3e52c9cf0adff3747c7683b69d371

<!-- end list -->

``` r
  digest::digest(mat, algo = 'md5')
```

    ## [1] "1fbbe9568738577a2f3e3dc42e6c75cf"

``` r
  digest::digest(metadata, algo = 'md5')
```

    ## [1] "542a40ccf8b14c51ffa45361c5d3aed9"

``` r
  digest::digest(sig_genes, algo = 'md5')
```

    ## [1] "fdc3e52c9cf0adff3747c7683b69d371"

Take a look at the contents of the data.

``` r
  # first 5 rows; first 5 columns
  mat[1:5,1:5]
```

    ##         SAM9103802 SAM9103803 SAM9103804 SAM9103805  SAM9103806
    ## A1BG     0.1288745  0.1637147 -0.1106011 -0.1113405  0.09776126
    ## A1CF     1.4491133  1.6378292  1.4676648  1.5119170  1.40215292
    ## A2M     15.0932787 14.8324464 14.7205315 14.6949978 14.70800150
    ## A2ML1    3.4826292  3.7443431  4.4786253  3.1842529  4.80886450
    ## A3GALT2  1.2259417  1.0704200  1.2452696  1.1006621  1.11011460

``` r
  # take a peek at the metadata
  head(metadata)
```

    ##            Pathotype CD3 CD20 CD68L CD68SL CD138
    ## SAM9103802  Lymphoid   2    3     1      3     3
    ## SAM9103803  Lymphoid   3    4     4      3     3
    ## SAM9103804   Myeloid   3    0     0      4     0
    ## SAM9103805  Lymphoid   2    2     3      1     1
    ## SAM9103806   Fibroid   0    0     0      1     0
    ## SAM9103807   Fibroid   0    0    NA      1     0

``` r
  # take a peek at the genes identified as statistically significant
  head(sig_genes)
```

    ## [1] "A2ML1"      "AADACL2"    "ABCA10"     "ABCA12"     "AC010646.3"
    ## [6] "ACACB"

``` r
  # dimensions of expression data and metadata, and length of sig_genes
  dim(mat)
```

    ## [1] 19279    87

``` r
  dim(metadata)
```

    ## [1] 87  6

``` r
  length(sig_genes)
```

    ## [1] 2772

``` r
  # verify integrity of metadata and expression matrix:
  # --check that both objects are aligned by name
  all(rownames(metadata) == colnames(mat))
```

    ## [1] TRUE

``` r
  # Subset the expression matrix for the statistically significant genes
  mat <- mat[sig_genes,]
```

# 4, generate the heatmap

#### \[a\] scale the data to Z-scores (by row)

This is quite standard when performing clustering and generating a
heatmap.

``` r
  heat <- t(scale(t(mat)))
```

#### \[b\] set colour scheme and choose breaks

``` r
  myCol <- colorRampPalette(c('dodgerblue', 'black', 'yellow'))(100)
  myBreaks <- seq(-3, 3, length.out = 100)
```

#### \[c\] create annotation: histo-pathotype and histology scores

First, we will just generate some colour mappings for the metadata
histology scores.

``` r
  # CD3
    cd3 <- metadata$CD3
    cd3 <- cd3[!is.na(cd3)] # remove missing values - we don't want to include these in the mapping
    pick.col <- brewer.pal(9, 'Greens')
    col.cd3 <- colorRampPalette(pick.col)(length(unique(cd3)))

  # CD20
    cd20 <- metadata$CD20
    cd20 <- cd20[!is.na(cd20)]
    pick.col <- brewer.pal(9, 'Blues')
    col.cd20 <- colorRampPalette(pick.col)(length(unique(cd20)))

  # CD68L
    cd68L <- metadata$CD68L
    cd68L <- cd68L[!is.na(cd68L)]
    pick.col <- brewer.pal(9, 'Reds')
    col.cd68L <- colorRampPalette(pick.col)(length(unique(cd68L)))

  # CD68SL
    cd68SL <- metadata$CD68SL
    cd68SL <- cd68L[!is.na(cd68L)]
    pick.col <- brewer.pal(9, 'Oranges')
    col.cd68SL <- colorRampPalette(pick.col)(length(unique(cd68SL)))

  # CD138
    cd138 <- metadata$CD138.max
    cd138 <- cd138[!is.na(cd138)]
    pick.col <- brewer.pal(9, 'Purples')
    col.cd138 <- colorRampPalette(pick.col)(length(unique(cd68SL)))
```

The use of *brewer.pal* and *colorRampPalette* above are to just
automatically produce hexidecimal colour codes that will be used later
for mapping between each histology score and each colour.

``` r
  unique(col.cd3)
```

    ## [1] "#F7FCF5" "#C7E9C0" "#74C476" "#238B45" "#00441B"

``` r
  unique(col.cd20)
```

    ## [1] "#F7FBFF" "#C6DBEF" "#6BAED6" "#2171B5" "#08306B"

``` r
  unique(col.cd68L)
```

    ## [1] "#FFF5F0" "#FCBBA1" "#FB6A4A" "#CB181D" "#67000D"

``` r
  unique(col.cd68SL)
```

    ## [1] "#FFF5EB" "#FDD0A2" "#FD8D3C" "#D94801" "#7F2704"

``` r
  unique(col.cd138)
```

    ## [1] "#FCFBFD" "#DADAEB" "#9E9AC8" "#6A51A3" "#3F007D"

Now let’s build the actual annotation object, i.e., the legend:

``` r
  # Create an initial data-frame of the annotation that we want to use
  # In this example, the 'ann' object turns out to be the exact same as 'metadata'
  ann <- data.frame(
    Pathotype = metadata$Pathotype,
    CD3 = metadata$CD3,
    CD20 = metadata$CD20,
    CD68L = metadata$CD68L,
    CD68SL = metadata$CD68SL,
    CD138 = metadata$CD138,
    stringsAsFactors = FALSE)

  # create the colour mapping
  colours <- list(
    Pathotype = c('Lymphoid' = 'blue', 'Myeloid' = 'red', 'Fibroid' = 'green3', 'Ungraded' = 'grey'),
    CD3 = c('0' = '#F7FCF5', '1' = '#C7E9C0', '2' = '#74C476', '3' = '#238B45', '4' = '#00441B'),
    CD20 = c('0' = '#F7FBFF', '1' = '#C6DBEF', '2' = '#6BAED6', '3' = '#2171B5', '4' = '#08306B'),
    CD68L = c('0' = '#FFF5F0', '1' = '#FCBBA1', '2' = '#FB6A4A', '3' = '#CB181D', '4' = '#67000D'),
    CD68SL = c('0' = '#FFF5EB', '1' = '#FDD0A2', '2' = '#FD8D3C', '3' = '#D94801', '4' = '#7F2704'),
    CD138 = c('0' = '#FCFBFD', '1' = '#DADAEB', '2' = '#9E9AC8', '3' = '#6A51A3', '4' = '#3F007D'))

  # now create the ComplexHeatmap annotation object
  # as most of these parameters are self-explanatory, comments will only appear where needed
  colAnn <- HeatmapAnnotation(
    df = ann,
    which = 'col', # 'col' (samples) or 'row' (gene) annotation?
    na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
    col = colours,
    annotation_height = 0.6,
    annotation_width = unit(1, 'cm'),
    gap = unit(1, 'mm'),
    annotation_legend_param = list(
      Pathotype = list(
        nrow = 4, # number of rows across which the legend will be arranged
        title = 'Pathotype',
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 12, fontface = 'bold'),
        labels_gp = gpar(fontsize = 12, fontface = 'bold')),
      CD3 = list(
        nrow = 5,
        title = 'CD3',
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 12, fontface = 'bold'),
        labels_gp = gpar(fontsize = 12, fontface = 'bold')),
      CD20 = list(
        nrow = 5,
        title = 'CD20',
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 12, fontface = 'bold'),
        labels_gp = gpar(fontsize = 12, fontface = 'bold')),
      CD68L = list(
        nrow = 5,
        title = 'CD68L',
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 12, fontface = 'bold'),
        labels_gp = gpar(fontsize = 12, fontface = 'bold')),
      CD68SL = list(
        nrow = 5,
        title = 'CD68SL',
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 12, fontface = 'bold'),
        labels_gp = gpar(fontsize = 12, fontface = 'bold')),
      CD138 = list(
        nrow = 5,
        title = 'CD138',
        title_position = 'topcenter',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 12, fontface = 'bold'),
        labels_gp = gpar(fontsize = 12, fontface = 'bold'))))
```

Believe me, there are many more parameters that can be configured than
those which I show here.

#### \[d\] create annotation: box-and-whisker plots

``` r
  boxplotCol <- HeatmapAnnotation(
    boxplot = anno_boxplot(
      heat,
      border = FALSE,
      gp = gpar(fill = '#CCCCCC'),
      pch = '.',
      size = unit(2, 'mm'),
      axis = TRUE,
      axis_param = list(
        gp = gpar(fontsize = 12),
        side = 'left')),
      annotation_width = unit(c(2.0), 'cm'),
      which = 'col')

  boxplotRow <- HeatmapAnnotation(
    boxplot = row_anno_boxplot(
      heat,
      border = FALSE,
      gp = gpar(fill = '#CCCCCC'),
      pch = '.',
      size = unit(2, 'mm'),
      axis = TRUE,
      axis_param = list(
        gp = gpar(fontsize = 12),
        side = 'top')),
      annotation_width = unit(c(2.0), 'cm'),
      which = 'row')
```

#### \[e\] create annotation: gene labels

Many heatmaps are produced from a large number of variables / genes,
which result in it being difficult to label each gene in the plot space.
Here, we can ‘step through’ the variables / genes and choose to only
label a select few.

The number of rows (genes) in our object is:

2772

In this code snippet, we ‘step through’ the rownames and only retain
every 40th successive label.

``` r
  genelabels <- rowAnnotation(
    Genes = anno_mark(
      at = seq(1, nrow(heat), 40),
      labels = rownames(heat)[seq(1, nrow(heat), 40)],
      labels_gp = gpar(fontsize = 10, fontface = 'bold'),
      padding = 0.75),
      width = unit(2.0, 'cm') +

      max_text_width(
        rownames(heat)[seq(1, nrow(heat), 40)],
        gp = gpar(fontsize = 10,  fontface = 'bold')))
```

#### \[f\] perform partitioning around medoids (PAM) to identify clusters in the data

Performing k-means or PAM on our data can help us to identify internal
‘structure’ in the data that may relate to biologically meaningful
pathways, as an example.

``` r
  pamClusters <- cluster::pam(heat, k = 4) # pre-select k = 4 centers
  pamClusters$clustering <- paste0('Cluster ', pamClusters$clustering)

  # fix order of the clusters to have 1 to 4, top to bottom
  pamClusters$clustering <- factor(pamClusters$clustering,
    levels = c('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4'))
```

#### \[g\] create the actual heatmap object

``` r
  hmap <- Heatmap(heat,

    # split the genes / rows according to the PAM clusters
      split = pamClusters$clustering,
      cluster_row_slices = FALSE,

    name = 'Gene\nZ-\nscore',

    col = colorRamp2(myBreaks, myCol),

    # parameters for the colour-bar that represents gradient of expression
      heatmap_legend_param = list(
        color_bar = 'continuous',
        legend_direction = 'vertical',
        legend_width = unit(8, 'cm'),
        legend_height = unit(5.0, 'cm'),
        title_position = 'topcenter',
        title_gp=gpar(fontsize = 12, fontface = 'bold'),
        labels_gp=gpar(fontsize = 12, fontface = 'bold')),

    # row (gene) parameters
      cluster_rows = TRUE,
      show_row_dend = TRUE,
      #row_title = 'Statistically significant genes',
      row_title_side = 'left',
      row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
      row_title_rot = 90,
      show_row_names = FALSE,
      row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
      row_names_side = 'left',
      row_dend_width = unit(25,'mm'),

    # column (sample) parameters
      cluster_columns = TRUE,
      show_column_dend = TRUE,
      column_title = '',
      column_title_side = 'bottom',
      column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
      column_title_rot = 0,
      show_column_names = FALSE,
      column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
      column_names_max_height = unit(10, 'cm'),
      column_dend_height = unit(25,'mm'),

    # cluster methods for rows and columns
      clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
      clustering_method_columns = 'ward.D2',
      clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
      clustering_method_rows = 'ward.D2',

    # specify top and bottom annotations
      top_annotation = colAnn,
      bottom_annotation = boxplotCol)
```

#### \[j\] draw the heatmap

``` r
  draw(hmap + genelabels,
    heatmap_legend_side = 'left',
    annotation_legend_side = 'right',
    row_sub_title_side = 'left')
```

![Molecular Portraits of Early Rheumatoid Arthritis Identify Clinical
and Treatment Response
Phenotypes](README_files/figure-gfm/clusterheatmap_fig1-1.png)

# 5, extra: change colour scheme, breaks, and do extra clustering on columns

Modifying colour and breaks can help to emphasise expression patterns in
the heatmap.

Columns (samples) can be clustered ‘on the fly’ via *column\_km*.

``` r
  myCol <- colorRampPalette(c('royalblue', 'white', 'red3'))(100)
  myBreaks <- seq(-1.5, 1.5, length.out = 100)

  hmap1 <- Heatmap(heat,

    name = 'Gene Z-score',

    col = colorRamp2(myBreaks, myCol),

    heatmap_legend_param = list(
      color_bar = 'continuous',
      legend_direction = 'horizontal',
      legend_width = unit(8, 'cm'),
      legend_height = unit(5.0, 'cm'),
      title_position = 'topcenter',
      title_gp=gpar(fontsize = 30, fontface = 'bold'),
      labels_gp=gpar(fontsize = 24, fontface = 'bold')),

    cluster_rows = TRUE,
    show_row_dend = TRUE,
    row_title = 'Statistically significant genes',
    row_title_side = 'left',
    row_title_gp = gpar(fontsize = 30,  fontface = 'bold'),
    row_title_rot = 90,
    show_row_names = FALSE,
    row_names_gp = gpar(fontsize = 11, fontface = 'bold'),
    row_names_side = 'left',
    row_dend_width = unit(25,'mm'),

    cluster_columns = TRUE,
    show_column_dend = TRUE,
    column_title = 'Samples',
    column_title_side = 'bottom',
    column_title_gp = gpar(fontsize = 30, fontface = 'bold'),
    column_title_rot = 0,
    show_column_names = FALSE,
    column_names_gp = gpar(fontsize = 8, fontface = 'bold'),
    column_names_max_height = unit(10, 'cm'),
    column_dend_height = unit(25,'mm'),

    clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
    clustering_method_columns = 'ward.D2',
    clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
    clustering_method_rows = 'ward.D2')


  myCol <- colorRampPalette(c('forestgreen', 'black', 'purple'))(100)
  myBreaks <- seq(-2, 2, length.out = 100)

  hmap2 <- Heatmap(heat,

    split = pamClusters$clustering,
    cluster_row_slices = FALSE,

    column_km = 6,

    name = 'Gene Z-score',
    col = colorRamp2(myBreaks, myCol),

    heatmap_legend_param = list(
      color_bar = 'continuous',
      legend_direction = 'horizontal',
      legend_width = unit(8, 'cm'),
      legend_height = unit(5.0, 'cm'),
      title_position = 'topcenter',
      title_gp=gpar(fontsize = 30, fontface = 'bold'),
      labels_gp=gpar(fontsize = 24, fontface = 'bold')),

    cluster_rows = TRUE,
    show_row_dend = FALSE,
    #row_title = 'Statistically significant genes',
    row_title_side = 'right',
    row_title_gp = gpar(fontsize = 30,  fontface = 'bold'),
    row_title_rot = 90,
    show_row_names = FALSE,
    row_names_gp = gpar(fontsize = 12, fontface = 'bold'),
    row_names_side = 'left',
    row_dend_width = unit(25,'mm'),

    cluster_columns = TRUE,
    show_column_dend = TRUE,
    column_title = 'Samples',
    column_title_side = 'bottom',
    column_title_gp = gpar(fontsize = 30, fontface = 'bold'),
    column_title_rot = 0,
    show_column_names = FALSE,
    column_names_gp = gpar(fontsize = 8, fontface = 'bold'),
    column_names_max_height = unit(10, 'cm'),
    column_dend_height = unit(25,'mm'),

    clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
    clustering_method_columns = 'ward.D2',
    clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
    clustering_method_rows = 'ward.D2')

  pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
      draw(hmap1,
        heatmap_legend_side = 'top',
        row_sub_title_side = 'left',
        newpage = FALSE)
      popViewport()

    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
      draw(hmap2,
        heatmap_legend_side = 'top',
        row_sub_title_side = 'right',
        newpage = FALSE)
      popViewport()
  popViewport()
```

![Molecular Portraits of Early Rheumatoid Arthritis Identify Clinical
and Treatment Response
Phenotypes](README_files/figure-gfm/clusterheatmap_fig2-1.png)

# 6, session info

``` r
sessionInfo()
```

    ## R version 3.6.3 (2020-02-29)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 16.04.6 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/atlas-base/atlas/libblas.so.3.0
    ## LAPACK: /usr/lib/atlas-base/atlas/liblapack.so.3.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=pt_BR.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=pt_BR.UTF-8    
    ##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=pt_BR.UTF-8   
    ##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] grid      stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ## [1] cluster_2.1.0        digest_0.6.25        circlize_0.4.8      
    ## [4] ComplexHeatmap_2.2.0 RColorBrewer_1.1-2   knitr_1.28          
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.4.6        magrittr_1.5        colorspace_1.4-1   
    ##  [4] clue_0.3-57         rjson_0.2.20        rlang_0.4.5        
    ##  [7] highr_0.8           stringr_1.4.0       tools_3.6.3        
    ## [10] parallel_3.6.3      xfun_0.13           png_0.1-7          
    ## [13] htmltools_0.4.0     yaml_2.2.1          BiocManager_1.30.10
    ## [16] GlobalOptions_0.1.1 shape_1.4.4         evaluate_0.14      
    ## [19] rmarkdown_2.1       stringi_1.4.6       compiler_3.6.3     
    ## [22] GetoptLong_0.1.8

# 7, references

R Core Team (2020)

Antoine Lucas et al. (2020)

Morgan (2019)

Maechler et al. (2019)

Lewis et al. (2019)

Gu, Eils, and Schlesner (2016)

Gu et al. (2014)

Neuwirth (2014)

<div id="refs" class="references">

<div id="ref-digest">

Antoine Lucas, Dirk Eddelbuettel with contributions by, Jarek Tuszynski,
Henrik Bengtsson, Simon Urbanek, Mario Frasca, Bryan Lewis, Murray
Stokely, et al. 2020. *Digest: Create Compact Hash Digests of R
Objects*. <https://CRAN.R-project.org/package=digest>.

</div>

<div id="ref-ComplexHeatmap">

Gu, Zuguang, Roland Eils, and Matthias Schlesner. 2016. “Complex
Heatmaps Reveal Patterns and Correlations in Multidimensional Genomic
Data.” *Bioinformatics*.

</div>

<div id="ref-circlize">

Gu, Zuguang, Lei Gu, Roland Eils, Matthias Schlesner, and Benedikt
Brors. 2014. “Circlize Implements and Enhances Circular Visualization in
R.” *Bioinformatics* 30 (19): 2811–2.

</div>

<div id="ref-E-MTAB-6141">

Lewis, Myles J, Michael R Barnes, Kevin Blighe, Katriona Goldmann,
Sharmila Rana, and et al. 2019. “Molecular Portraits of Early Rheumatoid
Arthritis Identify Clinical and Treatment Response Phenotypes.”
<https://pubmed.ncbi.nlm.nih.gov/31461658/.>

</div>

<div id="ref-cluster">

Maechler, Martin, Peter Rousseeuw, Anja Struyf, Mia Hubert, and Kurt
Hornik. 2019. *Cluster: Cluster Analysis Basics and Extensions*.

</div>

<div id="ref-BiocManager">

Morgan, Martin. 2019. *BiocManager: Access the Bioconductor Project
Package Repository*. <https://CRAN.R-project.org/package=BiocManager>.

</div>

<div id="ref-RColorBrewer">

Neuwirth, Erich. 2014. *RColorBrewer: ColorBrewer Palettes*.
<https://CRAN.R-project.org/package=RColorBrewer>.

</div>

<div id="ref-R">

R Core Team. 2020. *R: A Language and Environment for Statistical
Computing*. Vienna, Austria: R Foundation for Statistical Computing.
<https://www.R-project.org/>.

</div>

</div>
