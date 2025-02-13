import config
from pathlib import Path
import seq_analysis
import fast5_loader
import dorado_loader
import plot


def plot_all_rDNA_reads_in_dorado_bam(bam_path: str):
    with dorado_loader.BamToSingleReadReader(bam_path) as bam_reader:
        for single_reader in bam_reader:
            analyzer = seq_analysis.ReadAnalyzer(single_reader, config.RDNA_REF_HUMAN)
            if len(analyzer.loadeddata['sequence']) < 30000:
                continue
            split_analyzed_read = analyzer.split_alignment()
            if analyzer.check_rDNA_read_from_split_alignment():
                plot.signle_read_plot_structure(single_reader, split_analyzed_read, 'temp_figs/PSCA0047', met=True, pdf=False)


if __name__ == "__main__":
    bam_path = Path("test_files/201020_47/calls_2025-02-05_T09-56-59.bam")
    plot_all_rDNA_reads_in_dorado_bam(bam_path)