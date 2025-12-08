from typing import List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Bio import SeqIO
from alphagenome.visualization import plot_components
from alphagenome.models.dna_client import SUPPORTED_SEQUENCE_LENGTHS
from alphagenome.models.dna_client import OutputType

############################################################
# Output Parsing & Plotting
############################################################

def parse_output_types(output_type_names: List[str]):
    """Parse output type names from config into OutputType enums.

    Unsupported aliases are rejected with guidance.
    """

    output_types = []
    for name in output_type_names:
        name_upper = name.upper()
        if name_upper in {"SPLICE_JUNCTIONS", "SPLICE_SITES"}:
            raise ValueError(f"{name_upper} is not supported. Use SPLICE_SITE_USAGE instead.")
        if name_upper in OutputType.__members__:
            output_types.append(OutputType[name_upper])
        else:
            print(f"Warning: Unknown output type '{name}', skipping...")
    return output_types


def name_filter(names: List[str], substring: str) -> List[bool]:
    return [substring in n for n in names]


def filter_by_strand(track_data, strand: str):
    """Apply strand selection to a TrackData object."""
    s = str(strand).lower()
    mapping = {
        '+': track_data.filter_to_positive_strand,
        '-': track_data.filter_to_negative_strand,
        'nonnegative': track_data.filter_to_nonnegative_strand,
        'nonpositive': track_data.filter_to_nonpositive_strand,
        'stranded': track_data.filter_to_stranded,
        'unstranded': track_data.filter_to_unstranded,
    }
    return mapping.get(s, lambda: track_data)()


def plot_all_outputs(output, transcripts, requested_outputs, cfg, resize_width=2**15, target_interval=None, variant=None):
    """Plot requested outputs with optional track/strand filtering."""
    for output_type in requested_outputs:
        if variant:
            ref_data = output.reference.get(output_type)
            alt_data = output.alternate.get(output_type)
            if ref_data is None or alt_data is None:
                continue
            tdata = {'REF': ref_data, 'ALT': alt_data}
            colors = {'REF': 'dimgrey', 'ALT': 'red'}
        else:
            track_data = output.get(output_type)
            if track_data is None:
                continue
            tdata = track_data
            colors = None

        # Apply filters
        if cfg.get('track_filter'):
            if variant:
                tdata['REF'] = tdata['REF'].filter_tracks(name_filter(tdata['REF'].names, cfg['track_filter']))
                tdata['ALT'] = tdata['ALT'].filter_tracks(name_filter(tdata['ALT'].names, cfg['track_filter']))
            else:
                tdata = tdata.filter_tracks(name_filter(tdata.names, cfg['track_filter']))
        if 'strand' in cfg:
            if variant:
                tdata['REF'] = filter_by_strand(tdata['REF'], cfg['strand'])
                tdata['ALT'] = filter_by_strand(tdata['ALT'], cfg['strand'])
            else:
                tdata = filter_by_strand(tdata, cfg['strand'])

        components = [plot_components.TranscriptAnnotation(transcripts)]
        annotations=[plot_components.IntervalAnnotation(target_interval)]
        if variant:
            components.append(plot_components.OverlaidTracks(tdata=tdata, colors=colors))
            annotations  = [plot_components.VariantAnnotation([variant], alpha=0.6), plot_components.IntervalAnnotation(target_interval)]
        else:
            components.append(plot_components.Tracks(tdata, ylabel_template='{biosample_name} ({strand})\n{name}'))

        plot_components.plot(
            components,
            interval=tdata.interval.resize(resize_width) if not variant else tdata['REF'].interval.resize(resize_width),
            title=f'{output_type.name}: {cfg["gene_symbol"]}',
            annotations=annotations,
        )
    plt.tight_layout()
    plt.show()
    return annotations


############################################################
# ASO Window & Variants
############################################################

def find_exon_sequence(interval_start: int, exon_intervals):
    s, e = int(exon_intervals[0]), int(exon_intervals[1])
    return s - interval_start, e - interval_start


def enumerate_aso_variants(ref_sequence: str, start: int, end: int, aso_length: int):
    """Generate ASO-masked variants by sliding an 'N' window.
    Ensures the sliding window stays within [start, end) bounds.
    """
    n_block = 'N' * aso_length
    variants = [ref_sequence[:i] + n_block + ref_sequence[i + aso_length:] for i in range(start, end)]
    asos = [ref_sequence[i:i + aso_length] for i in range(start, end)]
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
    collapsed[mask] = ism_scores.sum(axis=1)
    return collapsed


############################################################
# Interval Sizing
############################################################

def get_optimal_resize_width(interval_width: int, config_resize: int | None = None) -> int:
    """Select an optimal resize width from supported lengths."""

    if config_resize is not None:
        return config_resize

    supported = sorted(SUPPORTED_SEQUENCE_LENGTHS.values())
    for length in supported:
        if length >= interval_width:
            return length

    max_len = supported[-1]
    print(f"⚠️ WARNING: Interval ({interval_width:,} bp) exceeds max supported ({max_len:,} bp). Using {max_len:,} bp; predictions will be truncated.")
    return max_len


############################################################
# Notebook Refactors: ASO scoring + BED export
############################################################

def score_asos(
    outputs,
    ref_output,
    requested_outputs,
    cfg,
    interval,
    variant_interval,
    asos,
    start: int,
    end: int,
    fasta_path: str,
):
    """Compute ASO scores and create SeqLogo plots.

    Returns (results_df, plots_list).
    """

    track_pattern = cfg.get('track_filter', '')

    plots = []
    results = pd.DataFrame(asos, columns=['ASO_sequence'])
    results['position'] = list(range(len(results)))

    ref_subseq = str(
        SeqIO.to_dict(SeqIO.parse(fasta_path, 'fasta'))[interval.chromosome]
        .seq[variant_interval.start:variant_interval.end]
    )
    mask = build_base_mask(ref_subseq)

    for output_type in requested_outputs:
        variant_array = []
        for i in range(len(outputs)):
            td = outputs[i].get(output_type)
            if td is None:
                continue
            if track_pattern:
                td = td.filter_tracks(name_filter(td.names, track_pattern))
            if 'strand' in cfg:
                td = filter_by_strand(td, cfg['strand'])
            variant_array.append(td.values.mean(axis=1))

        if not variant_array:
            print(f"Skipping {output_type.name}: no variant tracks")
            continue

        variant_array = np.stack(variant_array, axis=1)

        ref_td = ref_output.get(output_type)
        if ref_td is None:
            print(f"Skipping {output_type.name}: missing reference output")
            continue
        if cfg.get('track_filter'):
            ref_td = ref_td.filter_tracks(name_filter(ref_td.names, cfg['track_filter']))
        if 'strand' in cfg:
            ref_td = filter_by_strand(ref_td, cfg['strand'])

        ref_values = ref_td.values.mean(axis=1)[:, None]

        aso_scores = diff_mean(ref_values, variant_array, start, end)[:, None]
        collapsed = collapse_ism(aso_scores, mask)
        results[f'ASO_{output_type.name}'] = aso_scores.sum(1)

        ylabel = f"{output_type.name.lower()}"
        plots.append(
            plot_components.SeqLogo(
                scores=collapsed,
                scores_interval=variant_interval,
                ylabel=ylabel,
                max_width=cfg['flank'],
            )
        )

    return results, plots

def score_asos_and_export(
    outputs,
    ref_output,
    requested_outputs,
    cfg,
    interval,
    variant_interval,
    asos,
    samples_max,
    start: int,
    end: int,
    fasta_path: str,
    aso_length: int,
    results_dir: str,
    config_name: str,
):
    """Compute ASO scores, create SeqLogo plots, and export BED files.

    Returns (results_df, plots_list).
    """

    results, plots = score_asos(
        outputs=outputs,
        ref_output=ref_output,
        requested_outputs=requested_outputs,
        cfg=cfg,
        interval=interval,
        variant_interval=variant_interval,
        asos=asos,
        start=start,
        end=end,
        fasta_path=fasta_path,
    )

    for output_type in requested_outputs:
        col = f'ASO_{output_type.name}'
        if col not in results.columns:
            continue

        ref_td = ref_output.get(output_type)
        if ref_td is None:
            continue
        ref_score = ref_td.values[start:end].mean()

        top = results.sort_values(col, ascending=False).head(samples_max)
        top = top[top[col] > 0]
        bot = results.sort_values(col, ascending=True).head(samples_max)
        bot = bot[bot[col] < 0]

        def safe_div(a, b):
            return a / b if b not in (0, None) else 0.0

        def build_df(sorted_results, invert: bool):
            scores = sorted_results[col]
            # Normalize colors safely
            denom = ref_score if invert else (scores.max() if len(scores) else 1.0)
            denom = denom if denom not in (0, None) else 1.0
            norm_color = scores / denom
            vals = ((1 - norm_color) * 255 if not invert else (1 + norm_color) * 255).clip(1, 256)
            itemRgb = (
                [f'255,{int(v)},{int(v)}' for v in vals]
                if not invert else [f'{int(v)},{int(v)},255' for v in vals]
            )
            # Normalize score for label safely against ref_score
            norm_score = scores.apply(lambda s: safe_div(s, ref_score))
            # BED coordinates: position is 0-based relative to variant_interval.start
            chrom_starts = variant_interval.start + sorted_results['position']
            chrom_ends = chrom_starts + aso_length
            return pd.DataFrame({
                'chrom': [interval.chromosome] * len(sorted_results),
                'chromStart': chrom_starts,
                'chromEnd': chrom_ends,
                'name': [f"{'top' if not invert else 'bottom'}{i+1}({norm_score.iloc[i] * 100:.2f}%)" for i in range(len(sorted_results))],
                'score': scores,
                'strand': interval.strand,
                'thickStart': chrom_starts,
                'thickEnd': chrom_ends,
                'itemRgb': itemRgb,
            })

        bed_df = pd.concat([build_df(top, False), build_df(bot, True)])
        out_path = f"{results_dir}/{config_name}_ASO_{output_type.name}.bed"
        with open(out_path, 'w') as f:
            f.write(
                f'track name={config_name}_ASO_{output_type.name} description="ASO scores for {output_type.name}" visibility="pack" useScore=1\n'
            )
            bed_df.to_csv(f, sep='\t', index=False, header=False)

    return results, plots