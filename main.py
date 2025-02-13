import config
import seq_analysis
import fast5_loader
import plot


reader = fast5_loader.SingleFast5Reader("/mnt/data/nanopore_cas9_bc/221124_Y8585fN_CR_bc/workspace/0/00128cf9-6898-4859-a5dd-b94bb8b4d324.fast5")
analyzer = seq_analysis.ReadAnalyzer(reader, config.RDNA_REF_HUMAN)
split_analyzed_read = analyzer.split_alignment()
plot.signle_read_plot_structure(reader, split_analyzed_read, 'temp_figs/temp', met=True, pdf=False)