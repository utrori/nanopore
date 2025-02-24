import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import dorado_bam_io


def get_methlation_and_quality_from_bam(mapped_bam_path: str, dorado_bam_path: str, figname: str):
    with dorado_bam_io.BamAnalyzer(dorado_bam_path, mapped_bam_path) as analyzer:
        reads_data = analyzer.analyze()
        rid2quality = analyzer.quality_data_by_read_id
    all_met_scores = []
    quals = []
    coding_aves = []
    for read_id, read_data in reads_data.items():
        coding_mets = []
        ave_quality = np.mean(rid2quality[read_id])
        for alignment in read_data.alignments:
            met_list = alignment.methylation_data
            in_met_scores = []
            if 4000 < alignment.length < 16000:
                for i in met_list:
                    if alignment.query_alignment_start < i[0] < alignment.query_alignment_end:
                        in_met_scores.append(i[1])
                        coding_mets.append(i[1])
            if in_met_scores:
                all_met_scores.append(np.mean(in_met_scores))
        all_coding_met_ave = np.mean(coding_mets)
        quals.append(ave_quality)
        coding_aves.append(all_coding_met_ave)
    plt.scatter(quals, coding_aves, marker='.', alpha=0.4)
    plt.xlabel('Averaged phred quality score per read')
    plt.ylabel('Averaged methylation score per read')
    plt.savefig(figname, dpi=300)


if __name__ == '__main__':
    #get_methlation_and_quality_from_bam("test_files/HG02723_1/HG02723_1_dorado.bam", "test_files/HG02723_1/HG02723_1_mapped.bam", figname='temp_figs/HG02723_1_quality_vs_methylation.png')
    dorado_bam_path = "./test_files/dorado_output_PSCA0047/calls_2025-02-05_T09-56-59.bam"
    mapped_bam_path = Path("./test_files/201020_47_coding_mapped.bam")
    get_methlation_and_quality_from_bam(dorado_bam_path, mapped_bam_path, figname='temp_figs/201020_47_quality_vs_methylation.png')