import sys
import fanc
import logging
logging.basicConfig(level=logging.INFO)

input_file = sys.argv[1]
output_prefix = sys.argv[2]


def write_expected(input_file, output_prefix):
    hic = fanc.load(input_file, mode="a")
    intra_expected, expected_by_chromosome, inter_expected = hic.expected_values()
    bin_size = hic.bin_size
    distances = list(range(0, len(intra_expected) * bin_size, bin_size))

    output_file = output_prefix + "_expected_values_all.txt"

    with open(output_file, "w+") as out:
        out.write("distance\texpected\n")
        for i in range(0, len(intra_expected)):
            out.write(str(distances[i]) + "\t" + str(intra_expected[i]) + "\n")

    output_file = output_prefix + "_expected_values_per_chrom.txt"
    with open(output_file, "w+") as out:
        out.write("chrom\tdistance\texpected\n")
        for chrom in expected_by_chromosome.keys():
            expected = expected_by_chromosome[chrom]
            for i in range(0, len(expected)):
                out.write(chrom + "\t" + str(distances[i]) + "\t" + str(expected[i]) + "\n")


write_expected(input_file, output_prefix)
