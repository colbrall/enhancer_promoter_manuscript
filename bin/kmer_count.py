'''
kmer_count.py

for regions (FASTA file) counts all possible occurrences of sequences length k
and prints a bed file of their counts.

kmer_count.py k path/to/fasta
'''

import string
import sys
import itertools

def kmer_dict(k):
    d = {}
    out_string = "#chr\tstart\tend"
    for kmer in itertools.product('ACGT',repeat=k):
        d[string.join(kmer,'')] = 0
        out_string = out_string + "\t%s" % string.join(kmer,'')
    print out_string
    return d


def rev_comp(seq):
    new = ""
    for char in seq:
        if (char == "A"):
            new += "T"
        elif (char == "T"):
            new += "A"
        elif (char == "G"):
            new += "C"
        else:
            new += "G"
    return new[::-1]


def count(reg,dic,k):
    d = dic
    seq = ''
    ind = 0
    while ind < len(reg) and not reg[ind].startswith('>'):
		seq += reg[ind]
		ind += 1
    ind = 0
    while ind+k <= len(seq):
        kmer = ''
        for j in range(ind,ind+k):
            kmer += seq[j]
        try: d[str.upper(kmer)] += 1
        except: pass
        try: d[rev_comp(str.upper(kmer))] += 1
        except: pass
        ind += 1
    return d

def main():
    k = int(sys.argv[1])
    d = kmer_dict(k)
    reg_file = [line.strip() for line in open(sys.argv[2],'r')]
    for ind in range(len(reg_file)):
        if reg_file[ind].startswith('>'):
            d = count(reg_file[ind+1:],d,k)
            l = str.split(str.split(reg_file[ind].replace('>',''),':')[1],'-')
            out_string = "%s\t%s\t%s" % (str.split(reg_file[ind].replace('>',''),':')[0],l[0],l[1])
            for kmer in sorted(d.keys()):
                out_string += "\t%d" % d[kmer]
                d[kmer] = 0
            print out_string
        ind += 1

if __name__ == "__main__":
    main()
