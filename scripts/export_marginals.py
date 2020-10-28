import fanc
import sys
import logging
from fanc.tools.files import write_bed

logging.basicConfig(level=logging.INFO)

input_file = sys.argv[1]
output_file = sys.argv[2]

logging.info("Working on input file %s", str(input_file))


def export_marginals(hic_file, output_file):
    hic = fanc.load(hic_file)
    marginals = hic.marginals(masked=False, norm=False)
    regions = list(hic.regions())

    for pos, r in enumerate(regions):
        r.set_attribute("score", marginals[pos])

    write_bed(output_file, regions)


export_marginals(input_file, output_file)
