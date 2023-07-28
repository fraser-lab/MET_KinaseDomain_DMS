# MET_KD_fill_in_library.py
# G.Estevam


'''
- The The goal of this code is to produce an "fillin" codon oligo library for staturation mutagnesis of the MET kinase domain
- The script takes an input gene file, where the gene to be mutated is uppercase, and the flanking regions are lowercase
- Primers are created with a starting length of 40, with a NNK subsition per codon position
- If the primers have a Tm within the min/max Tm and have lengths within the min/max lengths, the oligos are saved
- Otherwise, nts are added or subtracted until it meets the requirements listed in the point above

input arg 2 MET_KD.txt
'''


'''
EXPECTED:

    5  gcaaaatactgtccacattgacctcagtgctctaAATCCAGAGCTGGTCCAGGCAGTGCAGCATGTAGTGATTGGGCCCAGTAGCCTGATT  3
            ||||||||||||||||||||||||||||
       <- 3 tatgacaggtgtaactggagtcacgagat 5 Reverse Primer

                                       5 NNKCCAGAGCTGGTCCAGGCAGTGCAGCATGTAGTGATT 3 -> Forward Primer
                                         ||||||||||||||||||||||||||||||||||||||
    3  cgttttatgacaggtgtaactggagtcacgagatTTAGGTCTCGACCAGGTCCGTCACGTCGTACATCACTAACCCGGGTCATCGGACTAA  5
'''



''' code start '''
import os
import sys
import math
import re
import csv
from Bio.SeqUtils import MeltingTemp as mt



''' global variables '''

sequence = open(sys.argv[1], "r")
seq = sequence.read()
protein_gene = ''.join([nt for nt in seq if nt.istitle()]) # uppercase sequence/ gene
a = re.search('[A-Z]', seq).start() # gene uppercase
gene = seq[a:] # gene + 3' lowercase flank
seq_length = len(seq)//3

max_primer_length = 51
min_primer_length = 30
max_Tm = 63
min_Tm = 61
gene_codons = len(gene)//3
protein_gene_codons = len(protein_gene)//3
forward_primers = []
reverse_primers = []

# Make a reverse lookup dictionary for codons and amino acids that only contains high-frequency codons
highfreqtable = {
	'L':'CTG',
	'I':'ATC',
	'M':'ATG',
	'F':'TTC',
	'V':'GTG',
	'P':'CCC',
	'T':'ACC',
	'A':'GCC',
	'Y':'TAC',
	'H':'CAC',
	'Q':'CAG',
	'N':'AAC',
	'K':'AAG',
	'D':'GAC',
	'E':'GAG',
	'C':'TGC',
	'*':'TGA',
	'W':'TGG',
	'S':'AGC',
	'R':'AGA',
	'G':'GGC',
}



''' functions '''

# complement of entire input sequence
def complement(seq):
    complement = {'a':'t', 'c':'g', 't':'a', 'g':'c','A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in seq])


# reverse complement of input sequence
# MissyEliot
# FlipItAndReverseIt
def reverse_complement(seq):
    complement = {'a':'t', 'c':'g', 't':'a', 'g':'c','A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in seq[::-1]])


# create primers for each codon substition at each position of the gene
# these are the forward primers, 5' - 3' direction, complementry to sense strand
def create_inverse_PCR_primer(seq, gene, gene_codons, protein_gene, NNK):
    for position in range(119, len(gene),3):
        codon = gene[position:position + 3]
        if codon.isupper():
            primer_length = 37
            primer_flank = gene[position+3 : position+primer_length]
            primer_name = "%d_FW_MET_KD_%s_mut" % ((position/3)+1,codon)
            primer = '%s%s' % (Stop, primer_flank)
            Tm = mt.Tm_GC(primer_flank) # Tm_GC module uses empirical formulas based on GC content
            if Tm < min_Tm:
                while Tm < min_Tm and len(primer) >= min_primer_length and len(primer) < max_primer_length:
                    primer_length += 1
                    primer_flank = gene[position+3 : position+primer_length]
                    Tm = mt.Tm_GC(primer_flank)
                    primer_name = "%d_FW_MET_KD_%s_mut" % ((position/3)+1,codon)
                    primer = '%s%s' % (Stop, primer_flank)
            if Tm > max_Tm:
                while Tm > max_Tm and len(primer) >= min_primer_length and len(primer) <= max_primer_length:
                    primer_length -= 1
                    primer_flank = gene[position+3 : position+primer_length]
                    Tm = mt.Tm_GC(primer_flank)
                    primer_name = "%d_FW_MET_KD_%s_mut" % ((position/3)+1,codon)
                    primer = '%s%s' % (Stop, primer_flank)
            print (primer_name, primer, Tm)
            fw_oligos = [primer_name, primer, Tm]
            forward_primers.append(fw_oligos)
    return forward_primers,primer_name, primer, Tm


# create reverse primers for each codon substition
# these are reverse complement primers, 5' - 3' direction, complementry to antisense strand
def create_reverse_primer (seq, gene, gene_codons, protein_gene):
    i = 1
    for position in range(0, len(seq),36):
        codon = seq[position:position + 3]
        if codon.isupper():
            codon = position
            reverse_primer_length = 30
            reverse_primer_name = "%d_RV_MET_KD_%s_mut" % (i,seq[codon:codon+3])
            reverse_primer = reverse_complement(seq[codon-reverse_primer_length : codon])
            rv_Tm = mt.Tm_GC(reverse_primer)
            if rv_Tm < min_Tm:
                while rv_Tm < min_Tm and len(reverse_primer) >= min_primer_length and len(reverse_primer) < max_primer_length:
                    reverse_primer_length += 1
                    reverse_primer_name = "%d_RV_MET_KD_%s_mut" % (i,seq[codon:codon+3])
                    reverse_primer = reverse_complement(seq[codon-reverse_primer_length : codon])
                    rv_Tm = mt.Tm_GC(reverse_primer)
            if rv_Tm > max_Tm:
                while rv_Tm > max_Tm and len(reverse_primer) >= min_primer_length and len(reverse_primer) <= max_primer_length:
                    reverse_primer_length -= 1
                    reverse_primer_name = "%d_RV_MET_KD_%s_mut" % (i,seq[codon:codon+3])
                    reverse_primer = reverse_complement(seq[codon-reverse_primer_length : codon])
                    rv_Tm = mt.Tm_GC(reverse_primer)
            i += 1
            print (reverse_primer_name, reverse_primer, rv_Tm)
            rv_oligos = [reverse_primer_name, reverse_primer, rv_Tm]
            reverse_primers.append(rv_oligos)
    return reverse_primers, reverse_primer_name, reverse_primer, rv_Tm

# create output file
def create_output_file (forward_primers,reverse_primers):

    fw_file = open('MET_KD_Y1194_fillin_library_fw.csv', 'w')
    writer = csv.writer(fw_file)
    writer.writerow(['Forward Primer Name ','Primer ', 'Tm'])
    writer.writerows(forward_primers)

    rv_file = open('MET_KD_Y1194_fillin_library_rv', 'w')
    writer2 = csv.writer(rv_file)
    writer2.writerow(['Reverse Primer Name ','Primer ', 'Tm'])
    writer2.writerows(reverse_primers)


''' call functions and print outputs '''

#print (gene)
#print (reverse_complement(seq))
#print (complement(seq))
create_inverse_PCR_primer(seq, gene, gene_codons, protein_gene, Stop)
create_reverse_primer (seq, gene, gene_codons, protein_gene)
print (seq)
print (len(gene))
create_output_file(forward_primers,reverse_primers)
