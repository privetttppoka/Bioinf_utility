from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqRecord import SeqRecord
from typing import Union
import os


class BiologicalSequence(ABC):
    """Abstract base class for biological sequences."""

    @abstractmethod
    def alphabet_chek(self):
        pass

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        return self.sequence[index]

    def __repr__(self):
        return f"{self.__class__.__name__}({self.sequence})"

    def __str__(self):
        return self.sequence


class NucleicAcidSequence(BiologicalSequence):
    """Base class for nucleic acid sequences (DNA and RNA)."""

    def __init__(self, sequence: str):
        self.sequence = sequence

    def alphabet_chek(self) -> bool:
        """Checks if the sequence contains only valid characters."""
        if not set(self.sequence).issubset(self.alphabet):
            raise ValueError("Invalid sequence")

    def complement(self):
        """Returns the complementary sequence. Must be implemented in derived classes."""
        raise NotImplementedError(
            "The complement method must be implemented in a subclass"
        )

    def reverse(self):
        """Returns the reverse of the sequence."""
        return self.__class__(self.sequence[::-1])

    def reverse_complement(self):
        """Returns the reverse complement of the sequence."""
        return self.__class__(self.reverse().complement())


class DNASequence(NucleicAcidSequence):
    """
    Class representing a DNA sequence.
    """

    def __init__(self, sequence: str, alphabet=set("ATGCatgc")):
        super().__init__(sequence)
        self.alphabet = alphabet

    def complement(self):
        """Returns the complementary DNA sequence."""
        complement_dict = {
            "A": "T",
            "T": "A",
            "G": "C",
            "C": "G",
            "a": "t",
            "t": "a",
            "g": "c",
            "c": "g",
        }
        return self.__class__("".join(complement_dict[base] for base in self.sequence))

    def transcribe(self):
        """Transcribes the DNA sequence into an RNA sequence."""
        return self.__class__(self.sequence.replace("T", "u").replace("t", "u"))


class RNASequence(NucleicAcidSequence):
    """
    Class representing an RNA sequence.
    """

    def __init__(self, sequence: str, alphabet=set("AUGCaugc")):
        super().__init__(sequence)
        self.alphabet = alphabet

    _complement_dict = {
        "A": "U",
        "U": "A",
        "G": "C",
        "C": "G",
        "a": "u",
        "u": "a",
        "g": "c",
        "c": "g",
    }


class AminoAcidSequence(BiologicalSequence):
    """
    Class representing an amino acid sequence.
    """

    def __init__(self, sequence: str):
        self.sequence = sequence

    def alphabet_chek(self):
        """Checks if the sequence contains only valid amino acid characters."""
        valid_aa_alphabet = set("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwya")
        if not set(self.sequence).issubset(valid_aa_alphabet):
            raise ValueError("Invalid amino acid sequence")

    def molecular_weight(self):
        """Calculates the molecular weight of the protein sequence."""
        aa_weights = {
            "A": 89.1,
            "C": 121.2,
            "D": 133.1,
            "E": 147.1,
            "F": 165.2,
            "G": 75.1,
            "H": 155.2,
            "I": 131.2,
            "K": 146.2,
            "L": 131.2,
            "M": 149.2,
            "N": 132.1,
            "P": 115.1,
            "Q": 146.2,
            "R": 174.2,
            "S": 105.1,
            "T": 119.1,
            "V": 117.1,
            "W": 204.2,
            "Y": 181.2,
        }
        self.molecular_weight = sum(
            aa_weights.get(aa.upper(), 0) for aa in self.sequence
        )
        return self.__class__(f"{self.molecular_weight:.2f}")


def filter_fastq(
    input_fastq: str,
    output_fastq: str = "",
    gc_bounds: Union[tuple[float, float], float] = (0, 100),
    length_bounds: Union[tuple[int, int], int] = (0, 2**32),
    quality_threshold: float = 0,
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

    output_dir = os.path.join(os.getcwd(), "filtered")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_fastq = os.path.join(output_dir, output_fastq)

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
            if (
                (length_bounds[0] <= sequence_length <= length_bounds[1])
                and (gc_bounds[0] <= gc_content <= gc_bounds[1])
                and (avg_quality >= quality_threshold)
            ):
                # Write the filtered record to the output file
                SeqIO.write(record, out_handle, "fastq")

