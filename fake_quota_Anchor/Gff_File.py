 #fix_gff and gff_file may be are correct

import numpy as np
import FastaFile
import re


class Transcript:
    cds_array = np.empty([0, 2], int)
    start = 0
    end = 0

    def __init__(self, name, strand, chr_name):
        self.name = name
        self.strand = strand
        self.chr_name = chr_name

    def add_cds(self, start, end):
        self.cds_array = np.append(self.cds_array, np.array([[start, end]]), axis=0)

    def update_order(self):
        if self.strand == "+":
            self.cds_array = np.sort(self.cds_array, axis=0)
            self.start = self.cds_array[0][0]
            self.end = self.cds_array[len(self.cds_array)-1][1]
        if self.strand == "-":
            self.cds_array = -np.sort(-self.cds_array, axis=0)
            self.start = self.cds_array[len(self.cds_array)-1][0]
            self.end = self.cds_array[0][1]

    def __lt__(self, other):
        if self.start < other.start:
            return True
        elif self.start == other.start and self.end < other.end:
            return True
        else:
            return False

    def __gt__(self, other):
        if self.start > other.start:
            return True
        elif self.start == other.start and self.end > other.end:
            return True
        else:
            return False

    def __eq__(self, other):
        if self.start == other.start and self.end == other.end:
            return True
        else:
            return False


class Gene:
    trans_array = np.empty([0, 1], Transcript)
    start = 0
    end = 0

    def __init__(self, name, strand, chr_name):
        self.name = name
        self.strand = strand
        self.chr_name = chr_name

    def add_trans(self, transcript):
        self.trans_array = np.append(self.trans_array, [[transcript]], axis=0)

    def update_order(self):
        self.trans_array = np.sort(self.trans_array, axis=0)
        self.start = self.trans_array[0][0].start
        for trans in self.trans_array:
            if trans[0].end > self.end:
                self.end = trans[0].end

    def __lt__(self, other):
        if self.start < other.start:
            return True
        elif self.start == other.start and self.end < other.end:
            return True
        else:
            return False

    def __gt__(self, other):
        if self.start > other.start:
            return True
        elif self.start == other.start and self.end > other.end:
            return True
        else:
            return False

    def __eq__(self, other):
        if self.start == other.start and self.end == other.end:
            return True
        else:
            return False


def read_gff(gff_file):     # "/home/xiaodong/Desktop/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff"
    chr_gene_trans_cds_dict = {}
    chr_trans_cds_dict = {}
    trans_gene_map = {}
    gene_array = np.empty([0, 1], Gene)
    gene_chr_dict = {}
    chr_gene_list = {}
    with open(gff_file) as f:
        for line in f:
            m = re.search("^(\S+)\t(\S+)\tCDS\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S)\t.+?;Parent=(.+?);", line)
            if m is not None:
                chr_name = m.group(1)
                start = int(m.group(3))
                end = int(m.group(4))
                strand = m.group(6)
                trans_name = m.group(8)
                if chr_name not in chr_trans_cds_dict:
                    chr_trans_cds_dict[chr_name] = {}
                if trans_name not in chr_trans_cds_dict[chr_name]:
                    chr_trans_cds_dict[chr_name][trans_name] = Transcript(trans_name, strand, chr_name)
                chr_trans_cds_dict[chr_name][trans_name].add_cds(start, end)   # trans attr cds_array
            else:
                m = re.search("^(\S+)\t(\S+)\tmRNA\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S)\tID=(.+?);Parent=(.+?);", line)
                if m is not None:
                    chr_name = m.group(1)
                    strand = m.group(6)
                    trans_name = m.group(8)
                    gene_name = m.group(9)
                    if gene_name not in trans_gene_map:
                        trans_gene_map[gene_name] = []
                    if gene_name not in gene_chr_dict:
                        gene_chr_dict[gene_name] = chr_name
                    trans_gene_map[gene_name].append(trans_name)
                    if chr_name not in chr_gene_trans_cds_dict:
                        chr_gene_trans_cds_dict[chr_name] = {}
                    if gene_name not in chr_gene_trans_cds_dict[chr_name]:
                        chr_gene_trans_cds_dict[chr_name][gene_name] = Gene(gene_name, strand, chr_name)
    for ch in chr_trans_cds_dict:
        for tran in chr_trans_cds_dict[ch]:
            chr_trans_cds_dict[ch][tran].update_order()  # dict includes attr of ordered cds_array, chr_name, start, end, strand.
        for ge in chr_gene_trans_cds_dict[ch]:
            for trans_name in trans_gene_map[ge]:
                transcript = chr_trans_cds_dict[ch][trans_name]
                chr_gene_trans_cds_dict[ch][ge].add_trans(transcript)
                chr_gene_trans_cds_dict[ch][ge].update_order()   # dict includes attr of ordered transcript, chr_name, start, end, strand.
            gene_array = np.append(gene_array, [[chr_gene_trans_cds_dict[ch][ge]]], axis=0)
            gene_array = np.sort(gene_array, axis=0)
            chr_gene_list[ch] = gene_array
    return chr_gene_trans_cds_dict, chr_gene_list, gene_chr_dict

# cds_seq is derived from cds sequence of transcript connecting
# trans_seq from the first cds start number to the last cds end


def update_cds_trans_seq(chr_gene_trans_cds_dict, chr_fasta_dict):
    for c in chr_gene_trans_cds_dict:
        for ge in chr_gene_trans_cds_dict[c]:
            i = 0
            strand = chr_gene_trans_cds_dict[c][ge].strand
            for trans in chr_gene_trans_cds_dict[c][ge].trans_array:
                start = trans[0].start
                end = trans[0].end
                chr_gene_trans_cds_dict[c][ge].trans_seq = FastaFile.sub_seq(chr_fasta_dict, c, start, end, strand)
                for cds in trans[0].cds_array:
                    st = cds[0]
                    en = cds[1]
                    cd = FastaFile.sub_seq(chr_fasta_dict, c, st, en, strand)
                    chr_gene_trans_cds_dict[c][ge].trans_array[i][0].cds_seq = ""
                    chr_gene_trans_cds_dict[c][ge].trans_array[i][0].cds_seq += cd
                i += 1

if __name__ == "__main__":
    chr_gene_dict, chr_gene_list, gene_chr_dict = read_gff("/home/songlab/Desktop/xd/genomic.gff")
    chr_names, chr_fasta = FastaFile.read_fasta_file("/home/songlab/Desktop/xd/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna")
    update_cds_trans_seq(chr_gene_dict, chr_fasta)
    print(chr_gene_dict["NC_050096.1"]["gene-LOC103644366"].trans_array[0][0].cds_seq)
