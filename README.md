### **filter_fastq and run_dna_rna_tools**

### Content:
* [For what?](#for-what)
* [fastq_module](#fastq_module)
* [run_dna_rna_tools](#run_dna_rna_tools)
* [List of procedures:](#list-of-procedures)
* [Developers](#developers)



## **For what?**

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

These functions are loaded from the bioinf_funcs module, which you can also use in your own programs via imports

## **Developers**
+ [Alina Nazarova](https://github.com/privetttppoka)
