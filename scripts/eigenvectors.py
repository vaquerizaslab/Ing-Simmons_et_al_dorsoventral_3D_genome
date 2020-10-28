import fanc
import os
import sys
import logging
logging.basicConfig(level=logging.INFO)

input_file = sys.argv[1]


def calc_2_eigenvectors(ab):
    eig1 = ab.eigenvector(eigenvector=0,
                          genome="/home/research/vaquerizas/store/genomes/insects/Dmel/6.07/fasta/dmel-all-chromosome-r6.07.fasta")
    eig2 = ab.eigenvector(eigenvector=1,
                          genome="/home/research/vaquerizas/store/genomes/insects/Dmel/6.07/fasta/dmel-all-chromosome-r6.07.fasta")

    return(eig1, eig2)


def write_bed_eig(cor_data, eig, file_name):
    with open(file_name, mode="w") as f:
        for i, region in enumerate(cor_data.regions):
            f.write("\t".join([str(region.chromosome), str(region.start),
                    str(region.end), str(float(eig[i]))]) + "\n")


def calc_and_write(hic_file):
    logging.info("working on %s", hic_file)

    eig1_out = hic_file.replace(".hic", "_eig1.bed")
    eig2_out = hic_file.replace(".hic", "_eig2.bed")
    if not os.path.exists(eig1_out):

        cor_file = hic_file.replace(".hic", ".cor")
        if os.path.exists(cor_file):
            logging.info("Correlation matrix %s exists, loading it", cor_file)
            ab = fanc.load(cor_file)
        else:
            logging.info("Calculating correlation matrix and saving to %s", cor_file)
            hic = fanc.load(hic_file)
            ab = fanc.ABCompartmentMatrix.from_hic(hic, file_name=cor_file)

        eig1, eig2 = calc_2_eigenvectors(ab)
        write_bed_eig(ab, eig1, file_name=eig1_out)
        write_bed_eig(ab, eig2, file_name=eig2_out)
    else:
        logging.info("output file exists; skipping!")


calc_and_write(input_file)
