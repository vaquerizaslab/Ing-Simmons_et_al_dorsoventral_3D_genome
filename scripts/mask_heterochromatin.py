import fanc
import shutil
from fanc.hic import HicEdgeFilter
import logging

logger = logging.getLogger(__name__)


class RegionFilter(HicEdgeFilter):
    def __init__(self, hic_object, region_based, mask=None):
        HicEdgeFilter.__init__(self, mask=mask)
        self.set_hic_object(hic)

        self._regions_to_mask = set()
        for filter_region in region_based.regions:
            for region in hic_object.regions(filter_region, lazy=True):
                self._regions_to_mask.add(region.ix)

        logger.info("Selected a total of {} regions to be masked".format(len(self._regions_to_mask)))

    def valid_edge(self, edge):
        """
        Check if an edge falls into a low-coverage region.
        """
        if edge.source in self._regions_to_mask:
            return False
        if edge.sink in self._regions_to_mask:
            return False
        return True


input_file = snakemake.input['hic']
output_file = snakemake.output[0]
mask_regions = snakemake.input['mask']

filter_regions = fanc.load(mask_regions)

shutil.copy(input_file, output_file)
hic = fanc.load(output_file, mode="a")

f = RegionFilter(hic, filter_regions)
hic.filter(f)
