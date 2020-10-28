import fanc
import logging
import sys
# import seaborn
from fanc.architecture.aggregate import loop_strength

logging.basicConfig(level=logging.INFO)

loops = sys.argv[1]
hic = sys.argv[2]
output_file = sys.argv[3]

loops = fanc.load(loops)  # must be BEDPE
hic = fanc.load(hic)

ls = loop_strength(hic, loops, pixels=5, keep_invalid=True)

with open(output_file, mode='w') as f:
    for i, l in enumerate(loops):
        f.write("\t".join(l.fields + [str(ls[i])]) + "\n")
