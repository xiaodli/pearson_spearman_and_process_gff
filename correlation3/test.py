import pandas as pd
from scipy.stats import pearsonr
list1 = pd.Series([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0225, 0, 0, 0])
list2 = pd.Series([0, 0, 0, 0, 0.022, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0.13, 0, 0, 0])
list3 = [1, 1, 1, 1]
list4 = [2, 2, 2, 2]
correlation_coefficient, _ = pearsonr(list1, list2)
print(correlation_coefficient)
#
# list1 = [(1, 2)]
# unzip = zip(*list1)

# import argparse
#
# # sub-command functions
#
#
# def foo(args):
#     print(args.x)
#     print(args.y)
#
#
# def bar(args):
#     print('((%s))' % args.z)
#
#
# # create the top-level parser
# parser = argparse.ArgumentParser()
# subparsers = parser.add_subparsers(required=True)
#
# # create the parser for the "foo" command
# parser_foo = subparsers.add_parser('foo')
# parser_foo.add_argument('-x', type=int, default=1)
# parser_foo.add_argument('-y', required=False)
# parser_foo.set_defaults(func=foo)
#
# # create the parser for the "bar" command
# parser_bar = subparsers.add_parser('bar')
# parser_bar.add_argument('z')
# parser_bar.set_defaults(func=bar)
#
# # parse the args and call whatever function was selected
# args = parser.parse_args()
# # if hasattr(args, 'func'):
# #     if '-x' in args.__dict__:
# args.func(args)
#     # elif args.corr_analysis:
#     #     corr_parser = subparsers.choices[args.corr_analysis]
#     #     corr_parser.print_help()
# # else:
# #     parser.print_help()
# print(args.func)
# print(args.__dict__)