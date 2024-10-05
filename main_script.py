from typing import Union
from additional_modules import bioinf_funcs


def filter_fastq(
    seqs: dict[str, tuple[str, str]],
    gc_bounds: Union[tuple[float, float], float] = (0, 100),
    length_bounds: Union[tuple[int, int], int] = (0, 2**32),
    quality_threshold: float = 0,
) -> dict[str, tuple[str, str]]:
    """The filter_fastq function filters fastq

    by CG-composition, length, and reading quality level.
    All these parameters can be changed by passing
    the necessary values to the function argument:
    gc_bounds(tuple),length_bounds(tuple), quality_threshold(float).
    """
    if isinstance(gc_bounds, float) or isinstance(gc_bounds, int):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, float) or isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)

    quality_summ = []
    final_dict = {}

    for name, (sequence, quality) in seqs.items():

        for char in quality:
            quality_summ.append(ord(char) - 33)
        cg_content = (sequence.count("C") + sequence.count("G")) * 100 / len(sequence)

        if (
            (length_bounds[0] <= len(sequence) <= length_bounds[1])
            and (gc_bounds[0] <= cg_content <= gc_bounds[1])
            and ((sum(quality_summ)) / len(quality_summ)) >= quality_threshold
        ):
            final_dict[name] = (sequence, quality)

    return final_dict


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
