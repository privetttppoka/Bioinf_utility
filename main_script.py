import os
from typing import Union
from additional_modules import bioinf_funcs


def fastq_module(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[tuple[float, float], float] = (0, 100),
    length_bounds: Union[tuple[int, int], int] = (0, 2**32),
    quality_threshold: float = 0,
) -> dict[str, tuple[str, str]]:
    """Filters sequences in a FASTQ file based on GC content, sequence length, and quality score.

    This function reads a FASTQ file, filters the sequences based on provided bounds for GC content,
    sequence length, and the quality threshold, then saves the filtered sequences into a new FASTQ file.
    """

    if isinstance(gc_bounds, float) or isinstance(gc_bounds, int):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, float) or isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)

    quality_summ = []
    folder_path = os.path.join(os.getcwd(), input_fastq)

    with open(folder_path, "r") as file:
        while True:
            name = file.readline().strip()
            if not name:
                break
            sequence = file.readline().strip()
            file.readline()
            quality = file.readline().strip()

            for char in quality:
                quality_summ.append(ord(char) - 33)
                cg_content = (
                    (sequence.count("C") + sequence.count("G")) * 100 / len(sequence)
                )

                if (
                    (length_bounds[0] <= len(sequence) <= length_bounds[1])
                    and (gc_bounds[0] <= cg_content <= gc_bounds[1])
                    and ((sum(quality_summ)) / len(quality_summ)) >= quality_threshold
                ):
                    output_dir = os.path.join(os.getcwd(), "filtered")

                    if not os.path.exists(output_dir):
                        os.makedirs(output_dir)

                    output_path = os.path.join(output_dir, output_fastq)

                    with open(output_path, "a") as file:
                        file.write(f"{name}\n{sequence}\n+\n{quality}\n")


def run_dna_rna_tools(*args: Union[str, list[str]]) -> Union[str, list[str]]:
    """The run_dna_rna_tools function can perform basic bioinformatic operations:

    transcribe, reverse_complement, reverse, complement.
    This function also checks for DNA or RNA sequences
    otherwise the sequence is not processed.
    Sequences are taken as input and the last argument
    is one of the operations to be performed on them:
    transcribe, reverse_complement, reverse, complement.
    """

    option = args[-1]
    seqs = list(args[:-1])

    seq_final = bioinf_funcs.check(seqs)

    if len(seq_final) != 0:

        if option == "reverse_complement":
            result = bioinf_funcs.reverse_complement(seq_final)
        elif option == "reverse":
            result = bioinf_funcs.reverse(seq_final)
        elif option == "transcribe":
            result = bioinf_funcs.transcribe(seq_final)
        elif option == "complement":
            result = bioinf_funcs.complement(seq_final)
        else:
            return "not an option"

    else:
        return "not dna or rna"

    if len(result) == 1:
        return result[0]
    return result
