import os
from Bio import SeqIO


def rev_comp(_sequence, _flag):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    if _flag == 1:
        # this is a text sequence and not a file
        return "".join(complement.get(nt, nt) for nt in reversed(_sequence))
    elif _flag == 2:
        # this is a file
        read_file = open(os.path.abspath(_sequence), "rU")
        out_rev_comp = open("out_rev_comp.fasta", "w")
        for seq_rec in SeqIO.parse(read_file,"fasta"):
            seq_id = seq_rec.id
            seq_seq = seq_rec.seq
            out_rev_comp.write(">"+seq_id+"\n")
            out_rev_comp.write("".join(complement.get(nt, nt) for nt in reversed(seq_seq))+"\n")