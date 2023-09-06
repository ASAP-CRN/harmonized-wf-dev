# Parse command-line arguments
library('argparse')
parser <- ArgumentParser(description='Plot QC metrics')
parser$add_argument('--working-dir', dest='working_dir', type='character', help='Working directory', default='/data/CARD_singlecell/harmony-rna/')
parser$add_argument('--script-dir', dest='script_dir', type='character', help='Directory containing workflow scripts', default='scripts')
parser$add_argument('--threads', dest='threads', type='integer', help='Number of threads to use for processing')
parser$add_argument('--metadata', dest='metadata', type='character', help='Metadata file output by gmm_doublet_calling')
parser$add_argument('--project-name', dest='project_name', type='character', help='Project name')
parser$add_argument('--output-violin-plots', dest='output_violin_plots', type='character', help='Name of output violin plot PDF')
parser$add_argument('--output-umis-genes-plot', dest='output_umis_genes_plot', type='character', help='Name of output umis vs. genes plot PDF')
args <- parser$parse_args()
script_dir <- args$script_dir

# Set working directory and load packages
setwd(args$working_dir)
source(paste0(script_dir, '/main/load_packages.r'))

# Set variables from args or snakemake parameters
threads <- if (is.null(args$threads)) snakemake@threads else args$threads
metadata <- if (is.null(args$metadata)) snakemake@input[['metadata']] else args$metadata
project_name <- if (is.null(args$project_name)) snakemake@params[['project_name']] else args$project_name
plot1_output_file <- if (is.null(args$output_violin_plots)) snakemake@output[['plot_1']] else args$output_violin_plots
plot2_output_file <- if (is.null(args$output_umis_genes_plot)) snakemake@output[['plot_2']] else args$output_umis_genes_plot

# Main
future::plan('multicore', workers=threads)
options(future.globals.maxSize=ngbs * 1000 * 1024^2)

m <- fread(metadata)

noise <- c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.rb', 'doublet_scores')

plot.list <- future.apply::future_lapply(noise, function(feature) {

    m %>%

        ggplot(aes(x=project_name, y=get(feature))) +

        geom_violin(fill='steelblue', color='black') + theme_bw() +

        theme(axis.text.x=element_blank(), axis.title.x=element_blank()) + ylab(feature)

})

p1 <- ggpubr::ggarrange(plotlist=plot.list, legend='none', align='hv', ncol=3, nrow=3)

ggsave(plot=p1, width=15, height=8, filename=plot1_output_file)
cat("Wrote QC plot 1:", plot1_output_file, "\n")

p2 <- m %>% ggplot(aes(x=get(noise[1]), y=get(noise[2]))) + geom_point(alpha=0.3) +
    theme_bw() + xlab('Number of UMIs') + ylab('Number of Genes')

ggsave(plot=p2, width=18, height=9, filename=plot2_output_file)
cat("Wrote QC plot 2:", plot2_output_file, "\n")