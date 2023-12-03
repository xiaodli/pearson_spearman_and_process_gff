import re


class Fasta:
    name = ""
    seq = []

    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

#            m = re.search("chromosome " + "(\\w)+", line)


def read_fasta_file(fastafile):
    chr_name = []
    chr_fasta_dict = {}
    string = ""
    with open(fastafile) as f:
        for line in f:
            m = re.search(">(\\S+)", line)
            if m is not None:
                name = m.group(1)
                chr_name.append(name)
                if (name is not None) and len(string) > 0:
                    chr_fasta_dict[chr_name[len(chr_name) - 2]] = Fasta(chr_name[len(chr_name) - 2], string)
                    string = ""
            else:
                s = re.sub("\\s", "", line)
                string += s.upper()
        chr_fasta_dict[chr_name[-1]] = Fasta(chr_name[-1], string)
    return chr_name, chr_fasta_dict


def reverse_complementary_seq(seq):
    seq = seq[::-1]
    sequence = ""
    for le in seq:
        if le == "A":
            le = "T"
        elif le == "G":
            le = "C"
        elif le == "C":
            le = "G"
        elif le == "T":
            le = "A"
        elif le == "S":
            le = "W"
        elif le == "M":
            le = "K"
        elif le == "K":
            le = "M"
        elif le == "Y":
            le = "R"
        elif le == "R":
            le = "Y"
        elif le == "W":
            le = "S"
        elif le == "B":
            le = "V"
        elif le == "V":
            le = "B"
        elif le == "H":
            le = "D"
        elif le == "D":
            le = "H"
        sequence += le
    return sequence


def sub_seq(chr_fasta_dict, name, start, end, strand):
    seq = chr_fasta_dict[name].seq
    if start > len(seq):
        return
    if end > len(seq):
        end = len(seq)
    sub = seq[start-1:end]
    if strand == "+":
        return sub
    else:
        sub = reverse_complementary_seq(sub)
        return sub


# Comma is essential for script running.  return error
# Is line equal to print(line) ?
#chr_name_list, character_fasta_dict = read_fasta_file("/home/xiaodong/Desktop/xd/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna")
#print(sub_seq(character_fasta_dict, "NC_050096.1", 500, 550, "+"))
#print(sub_seq(character_fasta_dict, "NC_050096.1", 500, 550, "-"))
