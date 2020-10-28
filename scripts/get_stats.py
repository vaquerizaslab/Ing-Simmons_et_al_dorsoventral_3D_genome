"""Get pairs filtering statistics for all samples in metadata.txt."""

import csv
import fanc
import logging
import os
import sys

logging.basicConfig(level=logging.INFO)

input_dir = sys.argv[1]
output_file = sys.argv[2]

logging.info("Working on input directory %s", input_dir)
logging.info("Will write to %s", output_file)

def stats(maskable, masked_table):
    """Get stats from pairs object and edges."""
    import tables as t
    statistics = maskable.mask_statistics(masked_table)
    # calculate total
    if isinstance(masked_table, t.Group):
        total = 0
        for table in masked_table:
            total += table._original_len()
    else:
        total = masked_table._original_len()
    return statistics, total


def write_stats(input_dir, output_file):
    """Write stats to file."""
    pairs_files = [fn for fn in os.listdir(input_dir) if fn.endswith("pairs")]
    for p in pairs_files:
        logging.info("Working on %s", p)
        pairs = fanc.load(os.path.join(input_dir, p))
        statistics, total = stats(pairs, pairs._edges)

        with open(output_file, "a+") as out:
            out = csv.writer(out)
            out.writerow([p, "total", total])
            for key, val in statistics.items():
                out.writerow([p, key, val])

write_stats(input_dir, output_file)
