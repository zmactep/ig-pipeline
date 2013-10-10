from Bio import SeqIO


def fasta_head(fasta_file, head_size):
    result = []
    for i, rec in enumerate(SeqIO.parse(fasta_file, "fasta")):
        if i < head_size:
            result.append(rec)
        else:
            break
    return result


def kmer_gen(seq, k):
    for i in range(len(seq) - k + 1):
        yield seq[i:i+k]
