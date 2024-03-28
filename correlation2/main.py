# import sys
import control_corr
import argparse
import process_to_R_boxplot_violin


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
    parser_sub1.add_argument("-e", "--expression", dest="expression", type=str, default="", help="gene expression matrix file")
    parser_sub1.add_argument("-co", "--corr", dest="corr", type=str, default="", help="corr result file")
    parser_sub1.add_argument("-m", "--method", dest="method", type=str, default="spearman", help="pearson or spearman")
    parser_sub1.add_argument("-b", "--bool", dest="bool_log", type=str, default="T", help="default log(tpm) if set F using tpm")
    parser_sub1.set_defaults(func=control_corr.row_corr)

    # subparser corr analysis process csv file to R input
    parser_sub2 = subparsers1.add_parser('process', help='process csv file to R boxplot_violin input')
    parser_sub2.add_argument("-i", "--input_file", dest="input", type=str, default="", help="six columns csv format file ")
    parser_sub2.add_argument("-o", "--output_file", dest="output", type=str, default="", help="boxplot_violin plot input file")
    parser_sub2.set_defaults(func=process_to_R_boxplot_violin.process)

    # block_corr analysis(third command)
    parser_sub3 = subparsers1.add_parser('block_corr', help='block gene corr analysis')
    parser_sub3.add_argument("-A", "--AnchorWave", dest="AnchorWave", type=str, default="", help="collinearity file by AnchorWave")
    parser_sub3.add_argument("-W", "--WGDI", dest="WGDI", type=str, default="", help="collinearity file by WGDI")
    parser_sub3.add_argument("-M", "--MCScanX", dest="MCScanX", type=str, default="", help="collinearity file by MCScanX")
    parser_sub3.add_argument("-J", "--jcvi", dest="jcvi", type=str, default="", help="collinearity file by jcvi")
    parser_sub3.add_argument("-e", "--expression", dest="expression", type=str, default="", help="gene expression matrix file")
    parser_sub3.add_argument("-co", "--corr", dest="corr", type=str, default="", help="corr result file")
    parser_sub3.add_argument("-m", "--method", dest="method", type=str, default="spearman", help="pearson or spearman")
    parser_sub3.add_argument("-c", "--column", dest="column", type=str, default="", help="e.g.typing R1-d1")
    parser_sub3.add_argument("-b", "--bool", dest="bool_log", type=str, default="T", help="default log(tpm) if set F using tpm")
    parser_sub3.add_argument("-B", "--Block_min", dest="Block_min", type=str, default="10",
                             help="using this value to restrict minimum block length for correlation analysis(default: 10) ")
    parser_sub3.add_argument("-r", "--is_random", dest="is_random", type=str, default="F",
                             help="random.shuffle(query_to_ref), to evaluate block is a functional region ? (T or F default：F)")
    parser_sub3.add_argument("-t", "--two_block", dest="two_block", type=str, default="F",
                             help="1+2 2+3 3+4 .... , to evaluate block is a functional region ? (T or F default：F)")
    parser_sub3.set_defaults(func=control_corr.block_total_corr)

    args = parser.parse_args()
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()
