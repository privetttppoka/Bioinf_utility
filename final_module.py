from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqRecord import SeqRecord
from typing import Union

def filter_fastq(
    input_fastq: str,
    output_fastq: str = "some_output",
    gc_bounds: Union[tuple[float, float], float] = (0, 100),
    length_bounds: Union[tuple[int, int], int] = (0, 2**32),
    quality_threshold: float = 0
) -> None:
    """Filters sequences in a FASTQ file based on GC content, sequence length, and quality score
    using Biopython.

    Args:
    - input_fastq (str): Input FASTQ file.
    - output_fastq (str): Output FASTQ file to store filtered sequences.
    - gc_bounds (tuple[float, float] or float): Bounds for GC content (default: (0, 100)).
    - length_bounds (tuple[int, int] or int): Bounds for sequence length (default: (0, 2**32)).
    - quality_threshold (float): Minimum average quality score (default: 0).
    """
    
    # If single values are provided for bounds, convert them to tuples
    if isinstance(gc_bounds, (float, int)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, (float, int)):
        length_bounds = (0, length_bounds)

    # Open output FASTQ file for writing filtered sequences
    with open(output_fastq, "w") as out_handle:
        # Iterate over the FASTQ file using Biopython's SeqIO
        for record in SeqIO.parse(input_fastq, "fastq"):
            sequence = record.seq
            quality = record.letter_annotations["phred_quality"]
            sequence_length = len(sequence)
            
            # Calculate GC content using Biopython's utility function
            gc_content = gc_fraction(sequence) * 100
            
            # Calculate average quality score
            avg_quality = sum(quality) / len(quality)
            
            # Apply the filtering conditions based on GC content, length, and quality
            if (length_bounds[0] <= sequence_length <= length_bounds[1]) and \
               (gc_bounds[0] <= gc_content <= gc_bounds[1]) and \
               (avg_quality >= quality_threshold):
                # Write the filtered record to the output file
                SeqIO.write(record, out_handle, "fastq")

