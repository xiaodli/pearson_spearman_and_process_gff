from argparse import ArgumentParser
import Gff_File
import sys


# query_gff  is path about query_gff
# subject_gff is path about ref_gff
def with_blast_protein_collinearity(query_gff, subject_gff, blast_result, outputfile):
    query_chromosome_to_gene_dict, query_chromosome_to_gene_list, query_gene_name_to_chr_dict = Gff_File.read_gff(query_gff)
    subject_chromosome_to_gene_dict, subject_chromosome_to_gene_list, subject_gene_name_to_chr_dict = Gff_File.read_gff(subject_gff)
    query_index_dict = dict()
    subject_index_dict = dict()
    for q in query_chromosome_to_gene_list:
        i = 1
        while i < len(query_chromosome_to_gene_list[q]):
            query_index_dict[query_chromosome_to_gene_list[q][i][0].name] = i
            i += 1
    for s in subject_chromosome_to_gene_list:
        i = 1
        while i < len(subject_chromosome_to_gene_list[s]):
            subject_index_dict[subject_chromosome_to_gene_list[s][i][0].name] = i
            i += 1
    match_pairs = set()
    out_ = open(outputfile, "w")
    with open(blast_result) as f:
        for line in f:
            elements = line.split()
            query_identity = elements[0]
            subject_identity = elements[1]
            similarity_percent = elements[2]
            alignment_length = int(elements[3])
            mismatch_length = elements[4]
            gap_opening = elements[5]
            query_start = elements[6]
            query_end = elements[7]
            subject_start = elements[8]
            send = elements[9]
            e_value = elements[10]
            bit_score = float(elements[11])
            match_one = query_identity+"_"+subject_identity
            if (match_one not in match_pairs) and (alignment_length > 250) and (bit_score > 250):
                match_pairs.add(match_one)
                out_.write(
                    query_identity + "\t" + str(query_gene_name_to_chr_dict[query_identity]) + "\t" +
                    str(query_index_dict[query_identity]) + "\t" +
                    str(query_chromosome_to_gene_dict[query_gene_name_to_chr_dict[query_identity]][query_identity].start) + "\t" +
                    str(query_chromosome_to_gene_dict[query_gene_name_to_chr_dict[query_identity]][query_identity].end) + "\t" +
                    query_chromosome_to_gene_dict[query_gene_name_to_chr_dict[query_identity]][query_identity].strand + "\t" +
                    subject_identity + "\t" + str(subject_gene_name_to_chr_dict[subject_identity]) + "\t" +
                    str(subject_index_dict[subject_identity]) + "\t" +
                    str(subject_chromosome_to_gene_dict[subject_gene_name_to_chr_dict[subject_identity]][subject_identity].start) + "\t" +
                    str(subject_chromosome_to_gene_dict[subject_gene_name_to_chr_dict[subject_identity]][subject_identity].end) + "\t" +
                    subject_chromosome_to_gene_dict[subject_gene_name_to_chr_dict[subject_identity]][subject_identity].strand + "\t" +
                    str(similarity_percent) + "\n")
                out_.close()


if __name__ == '__main__':
    parser = ArgumentParser(description="This is a format change about collinearity")
    parser.add_argument("-q", "--query",
                        help="query gff_file",
                        dest="query_gff",
                        default="")
    parser.add_argument("-r", "--ref",
                        help="ref_gff_file",
                        dest="ref_gff",
                        default="")
    parser.add_argument("-b", "--blast_result",
                        help="blast_result gff_file",
                        dest="blast_result",
                        default="")
    parser.add_argument("-o", "--output",
                        help="outputfile",
                        dest="output",
                        default="")
    args = parser.parse_args()
    if args.query_gff == "":
        print("Error occurs,input string is blank", file=sys.stderr)
        parser.print_help()
        sys.exit(1)
    if args.ref_gff == "":
        print("Error occurs,please specify --ref", file=sys.stderr)
        parser.print_help()
        sys.exit(1)
    if args.blast_result == "":
        print("Error occurs,please specify --blast_result", file=sys.stderr)
        parser.print_help()
        sys.exit(1)
    if args.output == "":
        print("Error occurs,input string is blank", file=sys.stderr)
        parser.print_help()
        sys.exit(1)
    query_gff_file = args.query_gff
    ref_gff_file = args.ref_gff
    blast_results = args.blast_result
    output = args.output    
    with_blast_protein_collinearity(query_gff_file, ref_gff_file, blast_results, output)


