### **FASTQ and DNA/RNA/Amino Acid Sequences Tool**

### Content:
* [For what?](#for-what)
* [DNASequence](#dnasequence)
* [RNASequence](#rnasequence)
* [AminoAcidSequence](#aminoacidsequence)
* [filter_fastq Function](#filter-fastq-function)
* [convert_multiline_fasta_to_oneline](#convert_multiline_fasta_to_oneline)
* [parse_blast_output](#parse_blast_output)
* [Developers](#developers)

## **For what?**
The *bioinf_modul_upgraded* contains several functions for working with **DNA, RNA, and Amino Acid sequences**, as well as **FASTQ and BLAST output files**.

## **DNASequence**
Class representing a DNA sequence.

### Possible Methods:
- `__len__` : Returns the length of the DNA sequence.
- `__str__` : Returns the DNA sequence as a string.
- `__repr__`: Prints the sequence in a readable format.
- **alphabet_check**: Checks if the sequence is a valid DNA sequence.
- **complement**: Returns the complementary DNA sequence.
- **transcribe**: Converts the DNA sequence into an RNA sequence.
- **reverse**: Returns the reverse of the DNA sequence.
- **reverse_complement**: Returns the reverse complement of the DNA sequence.

#### Usage example:
```python
from bioinf_modul_upgraded import DNASequence

dna = DNASequence('_your_sequence_here')
dna.complement()
dna.transcribe()
dna.reverse()
dna.reverse_complement()
```

## **RNASequence**
Class representing an RNA sequence.

### Possible Methods:
- `__len__` : Returns the length of the RNA sequence.
- `__str__` : Returns the RNA sequence as a string.
- `__repr__`: Prints the sequence in a readable format.
- **alphabet_check**: Checks if the sequence is a valid RNA sequence.
- **complement**: Returns the complementary RNA sequence.

#### Usage example:
```python
from bioinf_modul_upgraded import RNASequence

rna = RNASequence('_your_sequence_here')
rna.complement()
rna.reverse()
rna.reverse_complement()
```

## **AminoAcidSequence**
Class representing an amino acid sequence.

### Possible Methods:
- `__len__` : Returns the length of the amino acid sequence.
- `__str__` : Returns the amino acid sequence as a string.
- `__getitem__`: Returns a slice of the sequence.
- `__repr__`: Prints the sequence in a readable format.
- **alphabet_check**: Checks if the sequence contains only valid amino acid characters.
- **molecular_weight**: Computes the molecular weight of the sequence.

#### Usage example:
```python
from bioinf_modul_upgraded import AminoAcidSequence

AA = AminoAcidSequence('_your_sequence_here')
AA.molecular_weight()
AA.alphabet_check()
```

## **filter_fastq Function**
This function filters FASTQ data based on GC content, sequence length, and read quality threshold.

### Function arguments:
1. **`input_fastq: str`** - Path to the input FASTQ file.
2. **`output_fastq: str`** - Path to the output FASTQ file (default: saves in `filtered` folder).
3. **`gc_bounds: Union[tuple[float, float], float] = (0, 100)`** - GC content range.
4. **`length_bounds: Union[tuple[int, int], int] = (0, 2**32)`** - Sequence length range.
5. **`quality_threshold: float = 0`** - Minimum average quality score.

### Output:
Filtered sequences are saved in the `filtered` folder.

## **convert_multiline_fasta_to_oneline**
Converts a multiline FASTA file into a single-line format.

### Function arguments:
1. **`input_fasta: str`** - Path to the input FASTA file.
2. **`output_fasta: str`** (Optional) - Path to the output FASTA file. If not provided, the result is printed to the console.

### Usage example:
```python
from bioinf_modul_upgraded import convert_multiline_fasta_to_oneline

convert_multiline_fasta_to_oneline("multiline_fasta.fasta", "oneline_fasta.fasta")
```

## **parse_blast_output**
Extracts top BLAST hit descriptions for each query and writes them to an output file.

### Function arguments:
1. **`input_file: str`** - Path to the BLAST results file.
2. **`output_file: str`** - Path to the output file where the sorted list of top BLAST hits will be saved.

### Output:
The output file will contain one description per line, sorted alphabetically.

### Usage example:
```python
from bioinf_modul_upgraded import parse_blast_output

parse_blast_output("blast_results.txt", "sorted_hits.txt")
```

## **Developers**
+ [Alina Nazarova](https://github.com/privetttppoka)
