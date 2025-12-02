from typing import List
from pathlib import Path

import numpy as np


############################################################
# Output Parsing & Plotting
############################################################

def parse_output_types(output_type_names: List[str]):
    """Parse output type names from config into OutputType enums.

    Unsupported aliases are rejected with guidance.
    """
    from alphagenome.models.dna_client import OutputType

    output_types = []
    for name in output_type_names:
        name_upper = name.upper()
        if name_upper in {"SPLICE_JUNCTIONS", "SPLICE_SITES"}:
            raise ValueError(f"{name_upper} is not supported. Use SPLICE_SITE_USAGE instead.")
        try:
            output_types.append(OutputType[name_upper])
        except KeyError:
            print(f"Warning: Unknown output type '{name}', skipping...")
    return output_types


def name_filter(names: List[str], substring: str) -> List[bool]:
    return [substring in n for n in names]


def filter_by_strand(track_data, strand: str):
    """Apply strand selection to a TrackData object."""
    s = strand.lower()
    if s == '+':
        return track_data.filter_to_positive_strand()
    if s == '-':
        return track_data.filter_to_negative_strand()
    if s == 'nonnegative':
        return track_data.filter_to_nonnegative_strand()
    if s == 'nonpositive':
        return track_data.filter_to_nonpositive_strand()
    if s == 'stranded':
        return track_data.filter_to_stranded()
    if s == 'unstranded':
        return track_data.filter_to_unstranded()
    return track_data


def plot_all_outputs(output, transcripts, requested_outputs, cfg, resize_width=2**15, target_interval=None):
    """Plot requested outputs with optional track/strand filtering."""
    from alphagenome.visualization import plot_components
    import matplotlib.pyplot as plt

    for output_type in requested_outputs:
        track_data = output.get(output_type)
        if track_data is None:
            continue

        if cfg.get('track_filter'):
            track_data = track_data.filter_tracks(name_filter(track_data.names, cfg['track_filter']))
        if 'strand' in cfg:
            track_data = filter_by_strand(track_data, cfg['strand'])

        components = [plot_components.TranscriptAnnotation(transcripts),
                      plot_components.Tracks(track_data, ylabel_template='{biosample_name} ({strand})\n{name}')]

        plot_components.plot(
            components,
            interval=track_data.interval.resize(resize_width),
            title=f'{output_type.name}: {cfg["gene_symbol"]}',
            annotations=[plot_components.IntervalAnnotation(target_interval)]
        )
    plt.tight_layout()
    plt.show()


############################################################
# ASO Window & Variants
############################################################

def find_exon_sequence(interval_start: int, exon_intervals):
    s, e = int(exon_intervals[0]), int(exon_intervals[1])
    return s - interval_start, e - interval_start


def enumerate_aso_variants(ref_sequence: str, start: int, end: int, aso_length: int):
    """Generate ASO-masked variants by sliding an 'N' window."""
    variants, asos = [], []
    for i in range(start, end):
        masked_seq = ref_sequence[:i] + 'N' * aso_length + ref_sequence[i + aso_length:]
        variants.append(masked_seq)
        asos.append(ref_sequence[i:i + aso_length])
    return variants, asos


############################################################
# Scoring Utilities
############################################################

def diff_mean(ref: np.ndarray, alt: np.ndarray, start: int, end: int) -> np.ndarray:
    """Difference of means in [start:end]: alt - ref."""
    return alt[start:end].mean(0) - ref[start:end].mean(0)


def build_base_mask(sequence: str) -> np.ndarray:
    """One-hot mask for sequence bases (positions x 4)."""
    seq_np = np.array(list(sequence))
    return np.stack([seq_np == b for b in "ACGT"], -1)


def collapse_ism(ism_scores: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Collapse per-variant scores to position/base via mask."""
    collapsed = np.zeros(mask.shape, dtype=ism_scores.dtype)
    collapsed[mask] = ism_scores.sum(1)
    return collapsed


############################################################
# Interval Sizing
############################################################

def get_optimal_resize_width(interval_width: int, config_resize: int | None = None) -> int:
    """Select an optimal resize width from supported lengths."""
    from alphagenome.models.dna_client import SUPPORTED_SEQUENCE_LENGTHS

    if config_resize is not None:
        return config_resize

    supported = sorted(SUPPORTED_SEQUENCE_LENGTHS.values())
    for length in supported:
        if length >= interval_width:
            return length

    max_len = supported[-1]
    print(f"⚠️ WARNING: Interval ({interval_width:,} bp) exceeds max supported ({max_len:,} bp). Using {max_len:,} bp; predictions will be truncated.")
    return max_len