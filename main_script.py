from typing import Union
import bioinf_funcs

def filter_fastq(seqs: dict[str, tuple[str, str]],
                 gc_bounds: Union[tuple[float, float], float] = (0, 100),
                 length_bounds: Union[tuple[int, int], int] = (0, 2**32),
                 quality_threshold: float = 0) -> dict[str, tuple[str, str]]:
    
    quality_summ = []
    final_dict = {}
    
    for name, (sequence, quality) in seqs.items():
        
        for char in quality:
            quality_summ.append(ord(char) - 33)
            
        CG_content =  (sequence.count("C") + sequence.count("G")) * 100/ len(sequence)
        
        if (length_bounds[0]<= len(sequence) <= length_bounds[1]) and (gc_bounds[0]<=CG_content<=gc_bounds[1]) and ((sum(quality_summ))/len(quality_summ)) >= quality_threshold:
            final_dict[name] = (sequence, quality)  
            
    return final_dict