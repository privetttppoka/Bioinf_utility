list_dna_rna = "ATCGU"

dict_comp_dna = {
    "A": "T",
    "a": "t",
    "T": "A",
    "t": "a",
    "G": "C",
    "g": "c",
    "C": "G",
    "c": "g",
}
dict_comp_rna = {
    "A": "U",
    "a": "u",
    "U": "A",
    "u": "a",
    "G": "C",
    "g": "c",
    "C": "G",
    "c": "g",
}


def check(seqs):

    list_check = []

    for seq in seqs:
        sequ = seq.upper()
        if set(sequ) <= set(list_dna_rna) and not ("U" in sequ and "T" in sequ):
            list_check.append(seq)

    return list_check


def transcribe(seqs):

    return [seq.replace("T", "U").replace("t", "u") for seq in seqs]



def reverse(seqs):
    return [seq[::-1] for seq in seqs]


def complement(seqs):

    list_comp = []

    for seq in seqs:

        seq_new = ""

        dict_comp = dict_comp_rna if "U" in seq.upper() else dict_comp_dna
        seq_new = "".join([dict_comp[letter] for letter in seq])
        list_comp.append(seq_new)

    return list_comp


def reverse_complement(seqs):
    return reverse(complement(seqs))

