from Bio import SeqIO
import numpy as np
import plotly.graph_objs as go
import plotly.io as pio
import gzip
from plotly.subplots import make_subplots

pio.templates.default = "plotly"

# Generate read length data from a FASTQ file (returns tuple of two numpy arrays)
# If you have a list, you can skip this function and only convert the list to numpy array 

def calculate_n50(lengths: np.ndarray) -> int:
    """Calculate the N50 value from read lengths."""
    if len(lengths) == 0:
        return 0

    sorted_lengths = np.sort(lengths)[::-1]
    total = np.sum(sorted_lengths)
    cumulative = np.cumsum(sorted_lengths)
    n50_index = np.searchsorted(cumulative, total / 2)
    return sorted_lengths[n50_index]


def create_length_histogram(lengths: np.ndarray, sample_name: str, save_html: bool = True) -> go.Figure:
    """Create and return an interactive histogram plot for read lengths."""
    n50 = calculate_n50(lengths)
    max_length = np.max(lengths) if lengths.size > 0 else 0
    bins = max(round(int(max_length) / 500), 20)

    hist, bin_edges = np.histogram(lengths, bins=bins)

    # Plot the histogram
    fig = go.Figure()
    fig.add_trace(go.Bar(x=bin_edges[1:], y=hist,
                  marker_color='#4CB391', name="Read Lengths"))

    # Add N50 annotation
    if n50 > 0:
        fig.add_vline(x=n50, line_dash="dash", line_color="red")
        fig.add_annotation(
            text=f"N50: {n50} bp", x=n50, y=max(hist) * 0.9, showarrow=True, arrowhead=2
        )

    # Update layout
    fig.update_layout(
        title=f"{sample_name}: Read Length Distribution",
        xaxis_title="Read Length (bp)",
        yaxis_title="Frequency (read count)",
        template="plotly",
        margin=dict(t=50, b=50, l=50, r=50)
    )

    # Save to HTML if required
    if save_html:
        fig.write_html(f"{sample_name}_read_length_distribution.html")
        print(f"Saved: {sample_name}_read_length_distribution.html")

    return fig


def create_base_quality_histogram(read_qualities: np.array, sample_name: str, save_html: bool = True) -> go.Figure:
    """Create an interactive histogram for average base quality scores from a FASTQ file."""
    
    max_length = np.max(read_qualities) if read_qualities.size > 0 else 0
    bins = max(max(round(int(max_length) / 500), 10), 20)

    # Plot base quality histogram
    fig = go.Figure()
    hist, bin_edges = np.histogram(read_qualities, bins=bins)

    fig.add_trace(go.Bar(x=bin_edges[1:], y=hist,
                  marker_color='salmon', name="Base Quality"))

    # Update layout for base quality plot
    fig.update_layout(
        title=f"{sample_name}: Average Base Quality Distribution",
        xaxis_title="Base Quality Score (Phred)",
        yaxis_title="Frequency (read count)",
        template="plotly",
        margin=dict(t=50, b=50, l=50, r=50)
    )

    # Save to HTML if required
    if save_html:
        fig.write_html(f"{sample_name}_base_quality_distribution.html")
        print(f"Saved: {sample_name}_base_quality_distribution.html")

    return fig

def create_combined_histograms(read_lengths: list, read_qualities: list, sample_name: str, save_html: bool = True) -> None:
    """Create combined interactive histograms for read length and base quality and save them in one HTML file."""
    read_lengths = np.array(read_lengths)
    read_qualities = np.array(read_qualities)

    length_histogram = create_length_histogram(
        read_lengths, sample_name, save_html=False)
    base_quality_histogram = create_base_quality_histogram(
        read_qualities, sample_name, save_html=False)

    # Create subplots (2 columns)
    fig = make_subplots(
        rows=1, cols=2, subplot_titles=["Read Length Distribution", "Base Quality Distribution"],
        column_widths=[0.5, 0.5]
    )

    # Add traces for both histograms
    fig.add_trace(length_histogram['data'][0], row=1, col=1)
    fig.add_trace(base_quality_histogram['data'][0], row=1, col=2)

    # Add N50 annotation to the read length histogram
    n50 = calculate_n50(read_lengths)
    if n50 > 0:
        fig.add_vline(x=n50, row=1, col=1, line_dash="dash", line_color="purple")
        fig.add_annotation(
            text=f"N50: {n50} bp", x=n50, y=max(length_histogram['data'][0]['y']) * 0.9, showarrow=True, arrowhead=2,
            row=1, col=1
        )

    # Update layout
    fig.update_layout(
        title=f"{sample_name}: Read Length & Base Quality Distributions",
        showlegend=False,
        template="plotly",
        margin=dict(t=50, b=50, l=50, r=50)
    )

    # Save the combined plot to HTML
    if save_html:
        fig.write_html(f"{sample_name}_combined_histograms.html")
        print(f"Saved: {sample_name}_combined_histograms.html")


if __name__ == '__main__':
    path = '/home/gaox36/bioinfo_tools/on_target_LRS_pipeline/Transgene_Integration_LR_Pipeline/testdata/test_data/batch3_fq'
    fq_names = open(f'{path}/fq_filenames.txt').readlines()
    for fq in fq_names[:1]:
        fastq_path = f'{path}/{fq.strip()}'
        sample_name = fq.strip().split('.')[0]

        create_combined_histograms(fastq_path, sample_name)
        print(f'Finished {sample_name}')
