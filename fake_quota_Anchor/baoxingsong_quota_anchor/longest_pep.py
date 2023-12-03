# This method is effective for pep_names of one gene differing greatly. Not like as gene.1, gene.2, gene.3.
import numpy as np
from argparse import ArgumentParser
import fix_gff
import sys
import FastaFile


def longest_pep(gff_fil, fasta_fil, protein_file, longest_result):           # song deletes transcript array, leaving only one transcript through processing gff and fasta array.
    chr_gene_dict, chr_gene_list, gene_to_chr_dict = fix_gff.read_gff(gff_fil)
    chr_name_list, chr_fasta_dict = FastaFile.read_fasta_file(fasta_fil)     # Don't forget the position lacks two files.
    pep_name, pep_fasta_dict = FastaFile.read_fasta_file(protein_file)          # where is pep seq from ?
    fix_gff.update_cds_trans_seq(chr_gene_dict, chr_fasta_dict)
    for ch in chr_gene_list:
        for ge_na in chr_gene_dict[ch]:
            i = 0
            longest = 0
            for trans_one in chr_gene_dict[ch][ge_na].trans_array:   # tran_one is [tran] array.
                b = len(trans_one[0].cds_seq)
                if b > longest:
                    longest = b
                    if (len(chr_gene_dict[ch][ge_na].trans_array) >= 2) and (i >= 1):
                        np.delete(chr_gene_dict[ch][ge_na].trans_array, i-1, 0)
                        i -= 1
                if b < longest:
                    np.delete(chr_gene_dict[ch][ge_na].trans_array, i, 0)
                    i -= 1
                i += 1
    # trans_array exists one trans, processing over.

    with open(longest_result, "w")as f:
        for ch in chr_gene_list:
            for ge in chr_gene_dict[ch]:
                f.write(">"+ge+"\n")
                f.write(pep_fasta_dict[chr_gene_dict[ch][ge].trans_array[0][0].name].seq    +"\n")
    return chr_gene_dict


if __name__ == '__main__':
    parser = ArgumentParser(description="This is a get longest protein method but I think this exists problems")
    parser.add_argument("-g", "--gff_file",
                        help="gff_file",
                        dest="one_gff",
                        default="")
    parser.add_argument("-f", "--fasta_file",
                        help="fasta_file",
                        dest="fasta_file",
                        default="")
    parser.add_argument("-p", "--protein_file",
                        help="protein_file",
                        dest="protein_file",
                        default="")
    parser.add_argument("-o", "--output",
                        help="outputfile",
                        dest="output",
                        default="")
    args = parser.parse_args()
    if args.one_gff == "":
        print("Error occurs,input string is blank", file=sys.stderr)
        parser.print_help()
        sys.exit(1)
    if args.fasta_file == "":
        print("Error occurs,please specify --fasta_file", file=sys.stderr)
        parser.print_help()
        sys.exit(1)
    if args.protein_file == "":
        print("Error occurs,please specify --blast_result", file=sys.stderr)
        parser.print_help()
        sys.exit(1)
    if args.output == "":
        print("Error occurs,input string is blank", file=sys.stderr)
        parser.print_help()
        sys.exit(1)
    gff_file = args.one_gff
    fasta_file = args.fasta_file
    pep_file = args.protein_file
    output = args.output
    longest_pep(gff_file, fasta_file,  pep_file, output)

