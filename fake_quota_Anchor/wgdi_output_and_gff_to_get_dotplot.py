import sys

import pandas as pd
# /home/xiaodong/Desktop/AnchorWave/test/gold_danio_rerio
#   car_danio.collinearity  danio.wgdi.gff3   gold.wgdi.gff3

collinearity = pd.read_csv(sys.argv[1], sep="\t", comment="#", index_col=None)

