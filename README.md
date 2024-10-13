
### **filter_fastq and run_dna_rna_tools**

### Content:
* [For what?](#for-what)
* [parse_blast_output](#parse_blast_output)
* [convert_multiline_fasta_to_oneline](#convert_multiline_fasta_to_oneline)
* [fastq_module](#fastq_module)
* [run_dna_rna_tools](#run_dna_rna_tools)
* [List of procedures:](#list-of-procedures)
* [Developers](#developers)



## **For what?**
The *bio_files_processor.py*  contains two functions **parse_blast_output** and **convert_multiline_fasta_to_oneline** 

### parse_blast_output

The parse_blast_output function extracts the top BLAST hit descriptions for each query and writes them to an output file in a sorted order.

#### Function arguments:
1. 'input_file: str` Path to the BLAST results file.
2. `output_file: str` Path to the output file where the sorted list of top BLAST hits will be saved.

`Output:`
The output file will contain one description per line, sorted alphabetically.

#### Usage example 
```python
from bio_files_processor import parse_blast_output

# Set file paths
input_file = "example_blast_results.txt"
output_file = "sorted_blast_hits.txt"

# Parse BLAST output
parse_blast_output(input_file, output_file)
```

### convert_multiline_fasta_to_oneline

The convert_multiline_fasta_to_oneline function takes a FASTA file with sequences split across multiple lines and formats it so each sequence is presented on a single line.

#### Function arguments:
1. `input_fasta: str` Path to the input FASTA file with multiline sequences.
2. `output_fasta: str` (Optional) Path to the output FASTA file where single-line sequences will be saved. If not provided, the result will be printed to the console.
`Output:`
The output FASTA file will contain sequences in a single-line format.

#### Usage example
```python
from bio_files_processor import convert_multiline_fasta_to_oneline

# Set file paths
input_fasta = "multiline_fasta.fasta"
output_fasta = "oneline_fasta.fasta"

# Convert FASTA format
convert_multiline_fasta_to_oneline(input_fasta, output_fasta)
```

The main program *main_script.py* contains two functions **fastq_module** and **run_dna_rna_tools** that solve the bioinformatics problems described below: 

### fastq_module

This function **`fastq_module`** performs filtering of FASTQ data based on several parameters, including GC content, sequence length, and read quality threshold. Let's examine the operation of the function step by step:

#### Function arguments:
1. **`input_fastq: str`** is the name of the FASTQ file to be processed.
2. **`output_fastq: str`** - name of the output FASTQ file where the filtered data will be saved.
3. **`gc_bounds: union[tuple[float, float], float] = (0, 100)`** - bounds for GC content (GC-containing nucleotides: guanine and cytosine). You can pass either a tuple (lower and upper bounds) or a single value - then it will be interpreted as an upper bound and the lower bound will be 0.
4. **`length_bounds: union[tuple[tuple[int, int], int] = (0, 2**32)`** - bounds for the length of the sequence. Similarly, you can pass either a tuple or a single value - then it will be interpreted as the upper bound, and the lower bound will be 0.
5. **`quality_threshold: float = 0`** - minimum threshold for the average quality of the sequence.

   - If the sequence passed filtering, it is written to the output FASTQ file located in the `filtered` folder. If the folder does not exist, it is created automatically.

### Output:
This function filters the data from the FASTQ file based on three parameters:
1. **CG content**.
2. **Sequence length**.
3. **Average read quality**.

Filtered sequences are saved to a new FASTQ file in the `filtered` folder.


### run_dna_rna_tools

The run_dna_rna_tools takes as input an arbitrary number of arguments with DNA or RNA sequences (str), and the last argument is the name of the procedure to be executed. It then performs the specified action on all passed sequences and returns the result.

### List of procedures:

- transcribe - return the transcribed sequence
- reverse - return the expanded sequence
- complement - return the complementary sequence
- reverse_complement - return the reverse complementary sequence

### Usage example

```python
run_dna_rna_tools('ATG', 'transcribe') # 'AUG'
run_dna_rna_tools('ATG', 'reverse') # 'GTA'
run_dna_rna_tools('AtG', 'complement') # 'TaC'
run_dna_rna_tools('ATg', 'reverse_complement') # 'cAT'
run_dna_rna_tools('ATG', 'aT', 'reverse') # ['GTA', 'Ta']
```

These functions are loaded from the bioinf_funcs module, which you can also use in your own programs via imports

## **Developers**
+ [Alina Nazarova](https://github.com/privetttppoka)
