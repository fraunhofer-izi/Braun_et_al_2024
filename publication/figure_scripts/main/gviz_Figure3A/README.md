# Building several graphics and put them together as Figure 3A:

usage of R Gviz package (R studio 4.3.2, Gviz 1.47.1, GenomicRanges 1.54.1)

* ideograms: (01_Ideogram.R) chr3_KPNA4, chr20_ZGPAT, chr19_Polycomb-Associated_non_coding_RNAs
* full genes: (02_gviz_gene_insertion_site_function_mainTitle.R) KPNA4, ZGPAT, Polycomb-Associated_non_coding_RNAs
* zoomed in genes:(03_gviz_gene_insertion_site_function_zoom.R)  exons flanking integration Site
* as one plot: (04_gviz_local_integration_site_only_chim.R) custom genomic tracks depicting "intron - Carvykti - intron", coverage plots using datatracks(lines)
* adding lines (zoom effect) and coordiantes (exact IS) with Inkscape