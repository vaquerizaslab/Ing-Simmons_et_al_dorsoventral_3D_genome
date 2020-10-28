import fanc
import sys
import logging
logging.basicConfig(level=logging.INFO)

input_file = sys.argv[1]
threshold = int(sys.argv[2])
output_file = sys.argv[3]

logging.info("Working on input file {}, with threshold {}".format(input_file, str(threshold)))


def cis_local_ratio(hic, local_threshold=20000):
    """
    Calculate the cis/trans ratio for a Hic object.

    :param hic: :class:`~fanc,data.genomic.Hic` object
    :param local_threshold: Threshold, in base pairs, below which interactions
    will be considered 'local'.
    :return: tuple (ratio, local, far, factor)
    """
    cis_local = 0
    cis_far = 0
    regions_dict = hic.regions_dict
    for edge in hic.edges(lazy=True):
        if regions_dict[edge.source].chromosome == regions_dict[edge.sink].chromosome:
            distance = abs(regions_dict[edge.source].center - regions_dict[edge.sink].center)
            if distance < local_threshold:
                cis_local += edge.weight
            else:
                cis_far += edge.weight

    return cis_local / (cis_local + cis_far), cis_local, cis_far, 1.0


hic = fanc.load(input_file, mode='r')

r, cis, trans, f = cis_local_ratio(hic, local_threshold=threshold)

if output_file:
    with open(output_file, 'a') as o:
        o.write("file\tcis\ttrans\tratio\tfactor\tthreshold\n")
        o.write("{}\t{}\t{}\t{:.3f}\t{:.3f}\t{}\n".format(input_file, cis, trans, r, f, threshold))
print("{}".format(input_file))
print("\tcis: {}".format(cis))
print("\ttrans: {}".format(trans))
print("\tratio: {:.3f}".format(r))
print("\tfactor: {:.3f}".format(f))
hic.close()
