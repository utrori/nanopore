import matplotlib.pyplot as plt
import utilities
import config
import os
import numpy as np
import methylation_analysis
import fast5_loader
import dorado_bam_io
import dorado_loader

plt.rcParams['font.family']= 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'


def single_read_plot_structure_fast5(f5reader, analyzed_read, savedir, filename="", 
                                    met=True, adenine=False, pdf=False):
    """
    Plot read structure from a Fast5 file.
    
    Parameters:
    -----------
    f5reader : object
        Fast5 reader object
    analyzed_read : object
        Analyzed read object containing alignment information
    savedir : str
        Directory to save the plots
    filename : str, optional
        Custom filename for the plot (default: read_id)
    met : bool, optional
        Whether to include methylation data (default: True)
    adenine : bool, optional
        Whether to include adenine modification data (default: False)
    pdf : bool, optional
        Whether to save as PDF instead of PNG (default: False)
    """
    if not os.path.exists(savedir):
        os.makedirs(savedir, exist_ok=True)

    # Extract data from fast5
    loadeddata = f5reader.extract_read_data()
    read_id = loadeddata['read_id']
    seq = loadeddata['sequence']
    quality = loadeddata['quality_scores']
    methylation = loadeddata['methylation']
    
    # Determine species
    ref = analyzed_read.ref
    if 'human' in ref:
        species = 'human'
    elif 'mouse' in ref:
        species = 'mouse'
    elif 'fly' in ref:
        species = 'fly'
    else:
        species = 'unknown'
    
    split_align_res = analyzed_read.split_mapping_res
    
    # Use the common plotting function
    _plot_read_structure(read_id, seq, quality, methylation, split_align_res, 
                        species, savedir, filename, met, adenine, pdf)


def single_read_plot_structure_bam(read_id, reads_data, savedir, filename="", 
                                  met=True, adenine=False, pdf=False):
    """
    Plot read structure from pre-analyzed BAM data.
    
    Parameters:
    -----------
    read_id : str
        Read ID to plot
    reads_data : dict
        Dictionary of analyzed read data from BamAnalyzer.analyze()
    savedir : str
        Directory to save the plots
    filename : str, optional
        Custom filename for the plot (default: read_id)
    met : bool, optional
        Whether to include methylation data (default: True)
    adenine : bool, optional
        Whether to include adenine modification data (default: False)
    pdf : bool, optional
        Whether to save as PDF instead of PNG (default: False)
    """
    if not os.path.exists(savedir):
        os.makedirs(savedir, exist_ok=True)

    if read_id not in reads_data:
        print(f"Read {read_id} not found in the analyzed data")
        return
    
    read_data = reads_data[read_id]
    seq = read_data.sequence
    quality = read_data.get_quality_scores()
    
    if quality is None:
        print(f"No quality scores found for read {read_id}")
        return
    
    methylation = read_data.methylation_data
    
    # Determine species from reference name
    ref_name = read_data.alignments[0].reference_name if read_data.alignments else ""
    if 'human' in ref_name.lower():
        species = 'human'
    elif 'mouse' in ref_name.lower():
        species = 'mouse'
    elif 'fly' in ref_name.lower():
        species = 'fly'
    else:
        species = 'unknown'
    
    # Create split_align_res format for BAM data
    split_align_res = np.array([
        [
            alignment.aligned_segment.flag,
            '+' if not alignment.aligned_segment.is_reverse else '-',
            alignment.reference_start,
            alignment.aligned_segment.cigarstring,
            alignment.aligned_segment.get_tag('AS'),
            alignment.aligned_segment.query_length
        ]
        for alignment in read_data.alignments
    ])
    
    # Use the common plotting function
    _plot_read_structure(read_id, seq, quality, methylation, split_align_res, 
                        species, savedir, filename, met, adenine, pdf)


def single_read_plot_structure_from_bams(read_id, mapped_bam, dorado_bam, savedir, filename="", 
                                       met=True, adenine=False, pdf=False):
    """
    Plot read structure by analyzing BAM files.
    
    Parameters:
    -----------
    read_id : str
        Read ID to plot
    mapped_bam : str
        Path to mapped BAM file
    dorado_bam : str
        Path to Dorado BAM file with methylation data
    savedir : str
        Directory to save the plots
    filename : str, optional
        Custom filename for the plot (default: read_id)
    met : bool, optional
        Whether to include methylation data (default: True)
    adenine : bool, optional
        Whether to include adenine modification data (default: False)
    pdf : bool, optional
        Whether to save as PDF instead of PNG (default: False)
    """
    with dorado_bam_io.BamAnalyzer(mapped_bam, dorado_bam) as analyzer:
        reads_data = analyzer.analyze()
    
    single_read_plot_structure_bam(read_id, reads_data, savedir, filename, met, adenine, pdf)


def batch_plot_reads_from_bams(read_ids, mapped_bam, dorado_bam, savedir, met=True, adenine=False, pdf=False):
    """
    Plot structure for multiple reads from BAM files, efficiently reusing analysis.
    
    Parameters:
    -----------
    read_ids : list
        List of read IDs to plot
    mapped_bam : str
        Path to mapped BAM file
    dorado_bam : str
        Path to Dorado BAM file with methylation data
    savedir : str
        Directory to save the plots
    met : bool, optional
        Whether to include methylation data (default: True)
    adenine : bool, optional
        Whether to include adenine modification data (default: False)
    pdf : bool, optional
        Whether to save as PDF instead of PNG (default: False)
    """
    # Analyze BAM files once for all reads
    with dorado_bam_io.BamAnalyzer(mapped_bam, dorado_bam) as analyzer:
        reads_data = analyzer.analyze()
    
    # Process each read using the analyzed data
    for read_id in read_ids:
        if read_id in reads_data:
            single_read_plot_structure_bam(read_id, reads_data, savedir, 
                                         filename=read_id, met=met, adenine=adenine, pdf=pdf)
        else:
            print(f"Read {read_id} not found in BAM files")


def batch_plot_all_reads_from_bams(mapped_bam, dorado_bam, savedir, met=True, adenine=False, pdf=False, 
                                  max_reads=None, min_quality=None, skip_strand_switch=True):
    """
    Plot structure for all reads from BAM files efficiently.
    
    Parameters:
    -----------
    mapped_bam : str
        Path to mapped BAM file
    dorado_bam : str
        Path to Dorado BAM file with methylation data
    savedir : str
        Directory to save the plots
    met : bool, optional
        Whether to include methylation data (default: True)
    adenine : bool, optional
        Whether to include adenine modification data (default: False)
    pdf : bool, optional
        Whether to save as PDF instead of PNG (default: False)
    max_reads : int, optional
        Maximum number of reads to process (useful for large datasets)
    min_quality : float, optional
        Skip reads with average quality below this threshold
    skip_strand_switch : bool, optional
        Skip reads with strand switching (default: True)
    """
    # Analyze BAM files once for all reads
    with dorado_bam_io.BamAnalyzer(mapped_bam, dorado_bam) as analyzer:
        reads_data = analyzer.analyze()
    
    # Create output directory if it doesn't exist
    os.makedirs(savedir, exist_ok=True)
    
    print(f"Found {len(reads_data)} reads in BAM files")
    
    # Process reads
    processed_count = 0
    skipped_quality = 0
    skipped_strand_switch = 0
    
    for read_id, read_data in reads_data.items():
        # Check if we've reached the maximum number of reads
        if max_reads is not None and processed_count >= max_reads:
            print(f"Reached maximum number of reads to process ({max_reads})")
            break
        
        # Skip reads with low quality if specified
        if min_quality is not None:
            quality_scores = read_data.get_quality_scores()
            if quality_scores is None or np.mean(quality_scores) < min_quality:
                skipped_quality += 1
                continue
        
        # Skip reads with strand switching if specified
        if skip_strand_switch and read_data.has_strand_switch:
            skipped_strand_switch += 1
            continue
        
        # Plot the read
        single_read_plot_structure_bam(read_id, reads_data, savedir, 
                                     filename=read_id, met=met, adenine=adenine, pdf=pdf)
        processed_count += 1
        
        # Print progress every 10 reads
        if processed_count % 10 == 0:
            print(f"Processed {processed_count} reads...")
    
    # Print summary
    print(f"Plotting complete: {processed_count} reads plotted")
    if min_quality is not None:
        print(f"Skipped {skipped_quality} reads due to low quality (threshold: {min_quality})")
    if skip_strand_switch:
        print(f"Skipped {skipped_strand_switch} reads due to strand switching")
    print(f"All plots saved to {savedir}")


def _plot_read_structure(read_id, seq, quality, methylation, split_align_res, species, 
                        savedir, filename="", met=True, adenine=False, pdf=False):
    """
    Internal function to handle the actual plotting logic.
    
    Parameters:
    -----------
    read_id : str
        Read ID
    seq : str
        Read sequence
    quality : list
        Quality scores
    methylation : list
        Methylation data as list of (position, score) tuples
    split_align_res : np.ndarray
        Alignment results in the format used by utilities.get_read_structure
    species : str
        Species ('human', 'mouse', 'fly', or 'unknown')
    savedir : str
        Directory to save the plots
    filename : str, optional
        Custom filename for the plot (default: read_id)
    met : bool, optional
        Whether to include methylation data (default: True)
    adenine : bool, optional
        Whether to include adenine modification data (default: False)
    pdf : bool, optional
        Whether to save as PDF instead of PNG (default: False)
    """
    if not filename:
        filename = read_id
    
    fig = plt.figure()
    plt.subplots_adjust(left=0.2)
    ax = fig.add_subplot()
    met_window = 200
    
    # Read structure plot
    lc = utilities.get_read_structure(read_id, config.SPLIT_LENGTH, split_align_res, species=species)
    ax.add_collection(lc)
    
    # Quality plot
    bin1 = 200
    x = []
    y = []
    for n in range(0, len(seq) - bin1, bin1):
        x.append(n)
        y.append(np.mean(quality[n:n+bin1]) * 200 - 10000)
    ax.plot(x, y)
    
    # Methylation plot
    if met and methylation:
        x, y = methylation_analysis.methylation_binning(met_window, methylation, method='threshold')
        ax.bar(x, y, width=met_window, color='magenta', zorder=0, alpha=0.3)
    
    if adenine:
        # Handle adenine methylation if implemented
        pass
    
    # Set y-axis ticks and limits based on species
    if species == 'mouse' or species == 'human':
        ax.set_yticks((-10000, 0, 10000, 20000, 30000, 40000))
        ax.set_yticklabels(('UM', '0', '10k', '20k', '30k', '40k'), fontsize=20)
        ax.set_ylim([-12000, 47000])
    elif species == 'fly':
        ax.set_yticks((-2000, 0, 2000, 4000, 6000, 8000, 10000))
        ax.set_yticklabels(('UM', '0', '2k', '4k', '6k', '8k', '10k'), fontsize=20)
        ax.set_ylim([-2500, 14000])
    
    ax.set_title(f'Read: {read_id}, Mean quality: {np.mean(quality):.1f}', fontsize=15)
    ax.autoscale()
    
    # Ensure proper directory path format
    if savedir.endswith('/'):
        savedir = savedir[:-1]
    
    # Save the figure
    if pdf:
        plt.savefig(f'{savedir}/{filename}.pdf')
    else:
        plt.savefig(f'{savedir}/{filename}.png', dpi=300)
    plt.close()


# Keep this legacy function for backward compatibility
def signle_read_plot_structure(f5reader=None, analyzed_read=None, savedir=str, filename="", 
                              met=True, adenine=False, pdf=False,
                              dorado_bam=None, mapped_bam=None, read_id=None):
    """
    Legacy function that dispatches to the appropriate specialized function.
    This function is maintained for backward compatibility.
    
    Parameters:
    -----------
    f5reader : object, optional
        Fast5 reader object
    analyzed_read : object, optional
        Analyzed read object containing alignment information
    savedir : str
        Directory to save the plots
    filename : str, optional
        Custom filename for the plot (default: read_id)
    met : bool, optional
        Whether to include methylation data (default: True)
    adenine : bool, optional
        Whether to include adenine modification data (default: False)
    pdf : bool, optional
        Whether to save as PDF instead of PNG (default: False)
    dorado_bam : str, optional
        Path to Dorado BAM file with methylation data
    mapped_bam : str, optional
        Path to mapped BAM file with alignment data
    read_id : str, optional
        Read ID to plot when using BAM files
    """
    if f5reader and analyzed_read:
        # Fast5 approach
        single_read_plot_structure_fast5(f5reader, analyzed_read, savedir, filename, met, adenine, pdf)
    elif dorado_bam and mapped_bam and read_id:
        # BAM approach
        single_read_plot_structure_from_bams(read_id, mapped_bam, dorado_bam, savedir, filename, met, adenine, pdf)
    else:
        print("Error: Must provide either (f5reader and analyzed_read) or (dorado_bam, mapped_bam, and read_id)")