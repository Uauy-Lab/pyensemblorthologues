import pprint

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class EnsemblSequenceRegion:
    def __init__(self, data):
        self.data = data
        # pprint.pp(data)

    @property
    def seq(self):
        return self.data["seq"]

    @property
    def strand(self):
        return self.data["strand"]

    @property
    def start(self):
        return self.data["start"]

    @property
    def end(self):
        return self.data["end"]

    @property
    def species(self):
        return self.data["species"]

    @property
    def seq_region(self):
        return self.data["seq_region"]

    @property
    def description(self):
        return self.data["description"]

    @property
    def region(self):
        return f"{self.seq_region}:{self.start}-{self.end}"

    def __repr__(self) -> str:
        return f"<EnsembleSequenceRegion species:{self.species} region:{self.region} strand:{self.strand} >"

    def __len__(self):
        return abs(self.end - self.start)

    @property
    def record(self):
        record = SeqRecord(
            Seq(self.seq),
            id=f"{self.species}_{self.region}",
            name=self.species,
            description=f"{self.species} {self.region} {self.description}",
        )
        return record


class EnsemblPairwiseAlignment:
    def __init__(self, aln, base="triticum_aestivum"):
        self.aln = aln
        self.base_id = base
        self.__alignments = list(
            map(lambda a: EnsemblSequenceRegion(a), self.aln["alignments"])
        )

    @property
    def alignments(self):
        return self.__alignments

    @property
    def base(self):
        for a in self.alignments:
            if a.species == self.base_id:
                return a
        raise f"Unable to find base alignment {self}"

    @property
    def other(self):
        for a in self.alignments:
            if a.species != self.base_id:
                return a
        raise f"Unable to find other alignment {self}"

    def __len__(self):
        a = self.base
        return len(a)

    def __repr__(self):
        return (
            f"<EnsemblPairwiseAlignment len:{len(self)} alignments:{self.alignments} >"
        )


class EnsemblPairwiseAlignments:
    def __init__(self, response, base="triticum_aestivum"):
        self.response = response
        self.base = base
        self.alns = list(
            map(
                lambda aln: EnsemblPairwiseAlignment(aln, base=self.base), self.response
            )
        )

    def longest(self):
        ret = self.alns[0]
        print("Finding longest...")
        for aln in self.alns:
            pprint.pp(aln)
            if len(aln) > len(ret):
                ret = aln
        return ret
