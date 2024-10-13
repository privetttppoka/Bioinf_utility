import os


def parse_blast_output(input_file, output_file):
    """Parses a BLAST output file to extract and save the top hits for each query.

    This function reads the BLAST result file and extracts descriptions of the best alignments and
    from each significant alignment section saves the descriptions in a new file, sorted alphabetically.
    """

    input_dir = os.path.join(os.getcwd(), input_file)

    if not os.path.exists(input_dir):
        raise FileNotFoundError(f"Файл {input_file} не найден")

    best_hits = set()  # Используем множество для исключения дублирующихся записей

    with open(input_dir, "r") as infile:
        capture = False
        for line in infile:
            line = line.strip()
            if line.startswith("Sequences producing significant alignments:"):
                capture = True
            elif capture:
                if line == "" or line.startswith(">"):
                    capture = False
                else:
                    # Извлекаем описание из строки
                    description_1 = line.split("...")[0]
                    description = description_1.split("    ")[0]
                    best_hits.add(description)

    # Сортируем и записываем в выходной файл
    sorted_hits = sorted(best_hits)

    with open(output_file, "w") as outfile:
        for hit in sorted_hits:
            outfile.write(f"{hit}\n")


def convert_multiline_fasta_to_oneline(input_fasta, output_fasta):
    """Converts a multiline FASTA file into a single-line format for each sequence.

    This function reads a FASTA file in which sequences may be broken into multiple lines,
    concatenates the sequence lines for each header, and outputs the formatted FASTA file
    where each sequence appears in a single line after its header."""

    with open(input_fasta, "r") as file:

        sequences = {}
        header = None

        for line in file:
            line = line.strip()
            if line.startswith(">"):
                header = line
                sequences[header] = []
            else:
                sequences[header].append(line)

    formatted_sequences = []
    for header, seq_lines in sequences.items():
        formatted_sequences.append(f"{header}\n{''.join(seq_lines)}")

    if output_fasta:
        with open(output_fasta, "w") as outfile:
            outfile.write("\n".join(formatted_sequences))
    else:
        print("\n".join(formatted_sequences))
