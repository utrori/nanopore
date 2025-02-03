import config
import seq_analysis
import fast5_loader
import plot


reader = fast5_loader.SingleFast5Reader("/mnt/data/nanopore_cas9_bc/221124_Y8585fN_CR_bc/workspace/0/00128cf9-6898-4859-a5dd-b94bb8b4d324.fast5", methylation='cpg')
analyzer = seq_analysis.ReadAnalyzer(reader, config.RDNA_REF_HUMAN)
split_analyzed_read = analyzer.split_alignment()
plot.signle_read_plot_structure('temp_figs/temp', reader, split_analyzed_read.split_mapping_res, mouse='no', met='yes', pdf='no', adenine='no')