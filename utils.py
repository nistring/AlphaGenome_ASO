from typing import List
from pathlib import Path

import numpy as np



############################################################
# Track Filtering Helpers
############################################################

def name_filter(names: List[str], substring: str) -> List[bool]:
	return [substring in n for n in names]


def filter_by_strand(track_data, strand: str):
	"""Apply strand selection to a TrackData object.

	Supported values:
	'+'  -> positive strand only
	'-'  -> negative strand only
	'nonnegative' -> positive + unstranded
	'nonpositive' -> negative + unstranded
	'stranded' -> only stranded (+ or -)
	'unstranded' -> only unstranded '.'
	'all' -> no filtering
	"""
	strand = strand.lower()
	if strand == '+':
		return track_data.filter_to_positive_strand()
	if strand == '-':
		return track_data.filter_to_negative_strand()
	if strand == 'nonnegative':
		return track_data.filter_to_nonnegative_strand()
	if strand == 'nonpositive':
		return track_data.filter_to_nonpositive_strand()
	if strand == 'stranded':
		return track_data.filter_to_stranded()
	if strand == 'unstranded':
		return track_data.filter_to_unstranded()
	return track_data  # 'all'


def parse_output_types(output_type_names: List[str]):
	"""Parse output type names from config into OutputType enums.
	
	Args:
		output_type_names: List of output type names as strings (e.g., ['RNA_SEQ', 'SPLICE_SITES'])
	
	Returns:
		List of OutputType enum values
	
	Raises:
		ValueError: If SPLICE_JUNCTIONS is requested (not supported)
	"""
	from alphagenome.models.dna_client import OutputType
	
	output_types = []
	for name in output_type_names:
		name_upper = name.upper()
		if name_upper == 'SPLICE_JUNCTIONS' or name_upper == 'SPLICE_SITES':
			raise ValueError(f"{name_upper} is not supported. Use SPLICE_SITE_USAGE instead.")
		try:
			output_types.append(OutputType[name_upper])
		except KeyError:
			print(f"Warning: Unknown output type '{name}', skipping...")
	
	return output_types


def save_all_metadata(output, results_dir: str, gene_symbol: str):
	"""Save metadata for all available output types in the Output object.
	
	Args:
		output: AlphaGenome Output object containing predictions
		results_dir: Directory to save metadata files
		gene_symbol: Gene symbol for naming files
	"""
	from alphagenome.models.dna_client import OutputType
	
	gene_lower = gene_symbol.lower()
	results_path = Path(results_dir)
	results_path.mkdir(parents=True, exist_ok=True)
	
	output_type_names = {
		OutputType.ATAC: 'atac',
		OutputType.CAGE: 'cage',
		OutputType.DNASE: 'dnase',
		OutputType.RNA_SEQ: 'rna_seq',
		OutputType.CHIP_HISTONE: 'chip_histone',
		OutputType.CHIP_TF: 'chip_tf',
		OutputType.SPLICE_SITE_USAGE: 'splice_site_usage',
		OutputType.PROCAP: 'procap',
	}
	
	for output_type, name in output_type_names.items():
		track_data = output.get(output_type)
		if track_data is not None and hasattr(track_data, 'metadata'):
			filename = results_path / gene_lower / f"{name}_metadata.csv"
			track_data.metadata.to_csv(filename)
			print(f"Saved {output_type.name} metadata: {filename}")


def plot_all_outputs(output, transcripts, requested_outputs, cfg, resize_width=2**15):
	"""Plot all requested output types dynamically.
	
	Args:
		output: AlphaGenome Output object
		transcripts: Transcript annotations
		requested_outputs: List of OutputType enums
		cfg: Configuration dictionary
		resize_width: Width to resize intervals for plotting
	"""
	from alphagenome.visualization import plot_components
	import matplotlib.pyplot as plt
	

	for idx, output_type in enumerate(requested_outputs):
		track_data = output.get(output_type)
		if track_data is None:
			continue
		
		# Apply filtering based on config
		if 'track_filter' in cfg:
			track_data = track_data.filter_tracks(name_filter(track_data.names, cfg['track_filter']))
		if 'strand' in cfg:
			track_data = filter_by_strand(track_data, cfg['strand'])
		
		components = [plot_components.TranscriptAnnotation(transcripts)]
		components.append(plot_components.Tracks(track_data, ylabel_template='{biosample_name} ({strand})\n{name}'))
		
		plot_components.plot(
			components,
			interval=track_data.interval.resize(resize_width),
			title=f'{output_type.name}: {cfg["gene_symbol"]}',
		)
	plt.tight_layout()
	plt.show()


############################################################
# Variant Enumeration & Prediction
############################################################

def find_exon_sequence(ref_seq: str, exon_seq: str) -> tuple[int, int]:
	"""Find exon sequence in reference sequence and return coordinates.
	
	Args:
		ref_seq: Reference sequence to search
		exon_seq: Exon sequence to find
	
	Returns:
		Tuple of (start, end) coordinates relative to ref_seq
	
	Raises:
		ValueError: If sequence not found or found multiple times
	"""
	exon_seq_upper = exon_seq.upper()
	matches = []
	idx = 0
	while idx != -1:
		idx = ref_seq.find(exon_seq_upper, idx)
		if idx != -1:
			matches.append(idx)
			idx += 1
	
	if len(matches) == 0:
		raise ValueError(f"Exon sequence not found in reference sequence")
	if len(matches) > 1:
		raise ValueError(f"Exon sequence found {len(matches)} times in reference sequence at positions: {matches}")
	
	start = matches[0]
	return start, start + len(exon_seq_upper)


def enumerate_snv_variants(ref_sequence: str, start: int, end: int) -> List[str]:
	"""Generate all single-nucleotide variants for positions in [start, end).

	Coordinates are relative to ref_sequence index space.
	"""
	variants = []
	for i in range(start, end):
		for base in "ACGT":
			if base != ref_sequence[i]:
				variants.append(ref_sequence[:i] + base + ref_sequence[i + 1 :])
	return variants


def enumerate_aso_variants(ref_sequence: str, start: int, end: int, aso_length: int) -> List[str]:
	"""Generate ASO-masked variants by sliding a window of 'N's across positions.
	
	Args:
		ref_sequence: Reference sequence
		start: Start position (relative to ref_sequence)
		end: End position (relative to ref_sequence)
		aso_length: Length of ASO (number of bases to mask with 'N')
	
	Returns:
		List of sequences with ASO masking at different positions
	"""
	variants = []
	for i in range(start, end):
		masked_seq = ref_sequence[:i] + 'N' * aso_length + ref_sequence[i + aso_length:]
		variants.append(masked_seq)
	return variants

############################################################
# Scoring Utilities
############################################################

def diff_mean(ref: np.ndarray, alt: np.ndarray, start: int, end: int) -> np.ndarray:
	"""Compute difference of means across selected window."""
	return alt[start:end].mean(0) - ref[start:end].mean()


def build_base_mask(sequence: str) -> np.ndarray:
	"""Return one-hot mask for sequence bases (positions x 4)."""
	seq_np = np.array(list(sequence))
	return np.stack([seq_np == b for b in "ACGT"], -1)


def collapse_ism(ism_scores: np.ndarray, mask: np.ndarray) -> np.ndarray:
	"""Collapse ISM variant scores to per-position, per-base matrix using mask."""
	collapsed = np.zeros(mask.shape, dtype=ism_scores.dtype)
	collapsed[mask] = ism_scores.sum(1)
	return collapsed


def get_optimal_resize_width(interval_width: int, config_resize: int | None = None) -> int:
	"""Automatically select the optimal resize width for a given interval.
	
	Args:
		interval_width: Width of the gene interval
		config_resize: Optional resize width from config (takes precedence)
	
	Returns:
		Optimal resize width from supported sequence lengths
	"""
	from alphagenome.models.dna_client import SUPPORTED_SEQUENCE_LENGTHS
	
	# If config specifies a resize width, use it
	if config_resize is not None:
		return config_resize
	
	# Get sorted list of supported lengths
	supported_lengths = sorted(SUPPORTED_SEQUENCE_LENGTHS.values())
	
	# Find the smallest length that fits the interval
	for length in supported_lengths:
		if length >= interval_width:
			return length
	
	# If interval exceeds all supported lengths, use the largest and warn
	max_length = supported_lengths[-1]
	print(f"⚠️  WARNING: Interval width ({interval_width:,} bp) exceeds maximum supported length ({max_length:,} bp).")
	print(f"   Using maximum length {max_length:,} bp - predictions will be truncated.")
	return max_length