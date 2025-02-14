import numpy as np
import collections


def methylation_binning(window: int, methylation_data: list[(int, int)], method: str='average'):
    """Calculate average CpG methylation over a window. For plotting purposes."""
    threshold = 255 * 0.7
    summary_data = collections.defaultdict(list)
    x = []
    y = []
    for pos, score in methylation_data:
        summary_data[pos // window].append(score)
    max_pos = max(summary_data.keys())
    if method == 'average':
        for n in range(max_pos + 1):
            x.append(n * window)
            if summary_data.get(n, []):
                y.append(np.mean(summary_data.get(n, [])))
            else:
                y.append(0)
    elif method == 'threshold':
        for n in range(max_pos + 1):
            x.append(n * window)
            if summary_data.get(n, []):
                y.append(np.mean([1 if i > threshold else 0 for i in summary_data.get(n, [])])*10000)
            else:
                y.append(0)
    return x, y