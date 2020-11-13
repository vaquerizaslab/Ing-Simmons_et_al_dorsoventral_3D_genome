import os
import sys
import re
import fanc
from fanc.architecture.domains import InsulationScores
import logging
logging.basicConfig(level=logging.INFO)

input_file = sys.argv[1]


def write_insulation(hic_file):

    logging.info("working on %s", hic_file)
    hic = fanc.load(hic_file, mode='r')
    prefix = os.path.basename(hic_file).replace(".hic", "")

    res = int(re.findall('([0-9]+)kb', os.path.basename(hic_file))[0])
    logging.info("resolution detected as %s", str(res))

    window_sizes = [res * 1000 * w for w in [4, 6, 8, 10]]

    # calculate insulation index
    with InsulationScores.from_hic(hic, normalise=True, log=True, window_sizes=window_sizes,
                         file_name=os.path.join("data", "boundaries", prefix + "_micro-c.ii")) as ii:

        for window_size in window_sizes:
            w = window_size / res / 1000
            logging.info("Writing insulation index for window size %i", window_size)
            output_file = os.path.join("data", "boundaries", prefix + '_micro-c_{}.bw'.format(int(w)))
            ii.to_bigwig(output_file, window_size)


write_insulation(input_file)
