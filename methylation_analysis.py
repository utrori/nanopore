import numpy as np


def get_cpgs(seq):
    """Return a list of CpG sites."""
    cpg_sites = []
    upper_seq = seq.upper()
    for n in range(len(upper_seq) - 1):
        if upper_seq[n] == 'C':
            if upper_seq[n+1] == 'G':
                cpg_sites.append(n)
    return cpg_sites


def cpg_methylation_binning(window, seq, methylation_data, method='average'):
    """Calculate average CpG methylation over a window."""
    cpg_sites = get_cpgs(seq)
    threshold = 255 * 0.7
    x = []
    y = []
    if method == 'average':
        for n in range(0, len(seq) - window, window):
            x.append(n)
            in_cpgs = [i for i in cpg_sites if n < i < n + window]
            mod_scores = np.array(methylation_data)[in_cpgs]
            y.append(np.mean(mod_scores)*10000/255)
    elif method == 'threshold':
        for n in range(0, len(seq) - window, window):
            x.append(n)
            in_cpgs = [i for i in cpg_sites if n < i < n + window]
            mod_scores = np.array(methylation_data)[in_cpgs]
            y.append(np.mean([1 if i > threshold else 0 for i in mod_scores])*10000)
    return x, y