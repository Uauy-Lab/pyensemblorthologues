import os
import tempfile
from io import StringIO
from typing import List

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MafftCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class MSARegion:
    def __init__(self, regions=[], species="triticum_aestivum", reference=None):
        self.regions = regions
        self.species = species
        self.reference = reference

    # @property
    # def unaligned(self):
    #     ret = list()
    #     if len(self.regions) > 0:
    #         ret.append(self.regions[0].base.record)
    #         ret.extend(map(lambda a: a.other.record, self.regions))
    #     return ret

    def aligned(self):
        tmp_file = tempfile.NamedTemporaryFile(delete=False)
        print(tmp_file)
        print(tmp_file.name)
        tmp_file.close()
        SeqIO.write(self.regions, tmp_file.name, "fasta")
        mafft_cline = MafftCommandline(
            input=tmp_file.name
        )  # kwargs=["--localpair", "--quiet"
        # --ep 0 --genafpair
        # mafft_cline.localpair = True
        mafft_cline.genafpair = True
        # mafft_cline.ep = 0
        mafft_cline.quiet = True
        print(mafft_cline)
        stdout, stderr = mafft_cline()
        align = AlignIO.read(StringIO(stdout), "fasta")
        print(stderr)
        # tmp_file.delete()
        # os.remove(tmp_file.name)
        return align

    def add_sequence(self, sequence, id, name, description="", offset=0):
        sequence = ("-" * offset) + sequence
        seq = SeqRecord(sequence, id=id, name=name, description=description)
        self.regions.append(seq)
