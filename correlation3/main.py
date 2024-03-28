# import sys
import control_corr
import argparse

import specific_intersection_specific_all1_all2
import two_software_boxplot_violin
import two_software_venn

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='collinearity gene expression matrix correlation analysis, default: row correlation analysis')

    # subparser corr analysis
    subparsers1 = parser.add_subparsers(title='gene collinearity correlation analysis', dest='corr_analysis')

    # row_corr analysis(first command)
    parser_sub1 = subparsers1.add_parser('row_corr', help='single gene row corr analysis')
    parser_sub1.add_argument("-A", "--AnchorWave", dest="AnchorWave", type=str, default="", help="collinearity file by AnchorWave")
    parser_sub1.add_argument("-W", "--WGDI", dest="WGDI", type=str, default="", help="collinearity file by WGDI")
    parser_sub1.add_argument("-M", "--MCScanX", dest="MCScanX", type=str, default="", help="collinearity file by MCScanX")
    parser_sub1.add_argument("-J", "--jcvi", dest="jcvi", type=str, default="", help="collinearity file by jcvi")
    parser_sub1.add_argument("-t", "--tassel", dest="tassel", type=str, default="", help="maize tassel gene expression matrix file")
    parser_sub1.add_argument("-e", "--ear", dest="ear", type=str, default="", help=" maize ear gene expression matrix file")
    parser_sub1.add_argument("-s", "--sorghum", dest="sorghum", type=str, default="", help="sorghum gene expression matrix file")
    parser_sub1.add_argument("-co", "--corr", dest="corr", type=str, default="", help="corr result file")
    parser_sub1.add_argument("-m", "--method", dest="method", type=str, default="spearman", help="pearson or spearman")
    parser_sub1.add_argument("-b", "--bool", dest="bool_log", type=str, default="T", help="default log(tpm) if set F using tpm")
    parser_sub1.set_defaults(func=control_corr.row_corr)

    # transform format to get boxplot and violin by R
    parser_sub2 = subparsers1.add_parser('process_boxplot_violin', help='process csv file to R boxplot_violin input')
    parser_sub2.add_argument("-i", "--input_file", dest="input", type=str, default="", help="six columns csv format file ")
    parser_sub2.add_argument("-o1", "--output_file1", dest="output1", type=str, default="", help="boxplot_violin plot input file AnchorWave WGDI")
    parser_sub2.add_argument("-o2", "--output_file2", dest="output2", type=str, default="", help="boxplot_violin plot input file AnchorWave MCScanX")
    parser_sub2.add_argument("-o3", "--output_file3", dest="output3", type=str, default="", help="boxplot_violin plot input file AnchorWave jcvi")
    parser_sub2.set_defaults(func=two_software_boxplot_violin.process)

    # transform format to get venn by R
    parser_sub3 = subparsers1.add_parser('process_venn', help='process csv file to R venn input')
    parser_sub3.add_argument("-i", "--input_file", dest="input", type=str, default="", help="six columns csv format file ")
    parser_sub3.add_argument("-o1", "--output_file1", dest="output1", type=str, default="", help="boxplot_violin plot input file AnchorWave WGDI tassel")
    parser_sub3.add_argument("-o2", "--output_file2", dest="output2", type=str, default="", help="boxplot_violin plot input file AnchorWave WGDI ear")
    parser_sub3.add_argument("-o3", "--output_file3", dest="output3", type=str, default="", help="boxplot_violin plot input file AnchorWave MCScanX tassel")
    parser_sub3.add_argument("-o4", "--output_file4", dest="output4", type=str, default="", help="boxplot_violin plot input file AnchorWave MCScanX ear")
    parser_sub3.add_argument("-o5", "--output_file5", dest="output5", type=str, default="", help="boxplot_violin plot input file AnchorWave jcvi tassel")
    parser_sub3.add_argument("-o6", "--output_file6", dest="output6", type=str, default="", help="boxplot_violin plot input file AnchorWave jcvi ear")
    parser_sub3.set_defaults(func=two_software_venn.process)

    # transform format to get boxplot and violin
    parser_sub4 = subparsers1.add_parser('process_five', help='process csv file to R boxplot_violin five vs five column input')
    parser_sub4.add_argument("-i", "--input_file", dest="input", type=str, default="", help="six columns csv format file ")
    parser_sub4.add_argument("-o1", "--output_file1", dest="output1", type=str, default="", help="boxplot_violin plot input file AnchorWave WGDI tassel")
    parser_sub4.add_argument("-o2", "--output_file2", dest="output2", type=str, default="", help="boxplot_violin plot input file AnchorWave WGDI ear")
    parser_sub4.add_argument("-o3", "--output_file3", dest="output3", type=str, default="", help="boxplot_violin plot input file AnchorWave MCScanX tassel")
    parser_sub4.set_defaults(func=specific_intersection_specific_all1_all2.specific_all_inner)
    args = parser.parse_args()
    if hasattr(args, 'func'):
        # if "corr" in args.__dict__:
        args.func(args)
        # elif args.corr_analysis and ("corr" not in args.__dict__):
        # corr_parser = subparsers1.choices[args.corr_analysis]
        # corr_parser.print_help()
    else:
        parser.print_help()
