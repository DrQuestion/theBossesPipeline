# theBossesPipeline
A basic RNA-seq analysis pipeline, including differential expression analysis (dea) and functional enrichment analysis (fea).
It makes an intertwined use of bash and R, and it's callable from the Linux command line.
Made by two main bulky scripts: the first one, bigboss, it is developed to perform the computationally intensive tasks of the pipeline, the second one, theboss, to plot a peculiar version of the volcano plot and to perform fea of downregulated and upregulated genes. All the modular scripts in these two are written to be tunable under different options.

## bigboss
Performing the computationally intensive tasks, from the download of the fastqs to the dea, including quality check of the fastqs, alignment of the reads to the human genome and count of the features.
### downloader
### quality_checker
### aligner
### counter
### de_analyse.R

## theboss
Plots a volcano plot from the dea results and performs fea of GO.bp terms of the upregulated and the downregulated genes, plotting the results in a barplot representing the first ten terms and saving the whole tables in a tab separated format.
### plotte.R
### fe_analyse.R
