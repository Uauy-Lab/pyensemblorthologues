from typing import List

from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord


class MSARegion:
    def __init__(self, regions, species="triticum_aestivum"):
        self.regions = regions
        self.species = species

    @property
    def unaligned(self):
        ret = list()
        ret.append(self.regions[0].base.record)
        ret.extend(map(lambda a: a.other.record, self.regions))
        return ret

    def aligned(self):
        align = MultipleSeqAlignment(self.unaligned)
        return align
