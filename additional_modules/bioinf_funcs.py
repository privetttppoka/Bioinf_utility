list_DNA_RNA = "ATCGU"

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
        if set(sequ) <= set(list_DNA_RNA):
            list_check.append(seq)
        else:

    seq_check = []

    for seq in list_check:
        sequ = seq.upper()
        if (sequ.find("U") != -1) and (sequ.find("T") != -1):
            pass
        else:
            seq_check.append(seq)

    return seq_check




def transcribe(seqs):
    
    return [seq.replace("T", "U").replace("t", "u") for seq in seqs]


def reverse(seqs):

    list_rev = []

    for seq in seqs:
        list_rev.append(seq[::-1])

    return list_rev


def complement(seqs):

    list_comp = []
    for seq in seqs:
        seq_new = ""
        if (seq.find("u") == -1) and (seq.find("U") == -1):
            for letter in seq:
                seq_new = "".join([dict_comp_dna[letter] for letter in seq])
        else:
            for letter in seq:
                seq_new = "".join([dict_comp_rna[letter] for letter in seq])
        list_comp.append(seq_new)
    return list_comp


def reverse_complement(seqs):
    list_revcom = []
    seqs = complement(seqs)
    list_revcom = reverse(seqs)
    return list_revcom



