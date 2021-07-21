import pprint

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Mikado.parsers.GFF import GffLine

# from bitfield import Binary


class EnsemblSequenceRegion:
    def __init__(self, data):
        self.data = data
        # pprint.pp(data)

    @property
    def seq(self):
        return self.data["seq"]

    @property
    def strand(self):
        if self.data["strand"] > 0:
            return "+"
        elif self.data["strand"] < 0:
            return "-"
        return None

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

    @property
    def chrom(self):
        return self.seq_region

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

    @property
    def parent(self):
        return self.__parent

    @parent.setter
    def parent(self, parent):
        self.__parent = parent

    def __len__(self):
        a = self.base
        return len(a)

    def __repr__(self):
        return (
            f"<EnsemblPairwiseAlignment len:{len(self)} alignments:{self.alignments} >"
        )

    # def __cmp__(self, other):
    #     base  = self.base
    #     other_base = other.base

    #     non_ref = self.other
    #     non_ref_other = other.other

    #     if other_base.chrom != base.chrom:
    #         return base.chrom.__cmp__(other_base.chrom)
    #     if other_base.start != base.start:
    #         return base.start.__cmp__(other_base.start)
    #     if other_base.end != base.end:
    #         return base.end.__cmp__(other_base.end)
    #     if non_ref.species != non_ref_other.species:
    #         return non_ref.species.__cmp__(non_ref_other.species)
    #     return non_ref.__repr__().__cmp__(non_ref_other.__repr__())

    def __eq__(self, other):
        base = self.base
        other_base = other.base
        non_ref = self.other
        non_ref_other = other.other
        return (base.chrom, base.start, base.end, non_ref.species) == (
            other_base.chrom,
            other_base.start,
            other_base.end,
            non_ref_other.species,
        )

    def __ne__(self, other):
        base = self.base
        other_base = other.base
        non_ref = self.other
        non_ref_other = other.other
        return (base.chrom, base.start, base.end, non_ref.species) != (
            other_base.chrom,
            other_base.start,
            other_base.end,
            non_ref_other.species,
        )

    def __lt__(self, other):
        base = self.base
        other_base = other.base
        non_ref = self.other
        non_ref_other = other.other
        return (base.chrom, base.start, base.end, non_ref.species) < (
            other_base.chrom,
            other_base.start,
            other_base.end,
            non_ref_other.species,
        )

    def __le__(self, other):
        base = self.base
        other_base = other.base
        non_ref = self.other
        non_ref_other = other.other
        return (base.chrom, base.start, base.end, non_ref.species) <= (
            other_base.chrom,
            other_base.start,
            other_base.end,
            non_ref_other.species,
        )

    def __gt__(self, other):
        base = self.base
        other_base = other.base
        non_ref = self.other
        non_ref_other = other.other
        return (base.chrom, base.start, base.end, non_ref.species) > (
            other_base.chrom,
            other_base.start,
            other_base.end,
            non_ref_other.species,
        )

    def __ge__(self, other):
        base = self.base
        other_base = other.base
        non_ref = self.other
        non_ref_other = other.other
        return (base.chrom, base.start, base.end, non_ref.species) >= (
            other_base.chrom,
            other_base.start,
            other_base.end,
            non_ref_other.species,
        )

    def gff(self, seq=False):
        base = self.base
        other = self.other
        # print(other)
        gff = GffLine("")
        gff.attributes["chrom"] = base.seq_region
        gff.attributes["start"] = base.start
        gff.attributes["end"] = base.end
        gff.attributes["strand"] = base.strand
        gff.attributes["source"] = "Ensembl_compara"
        gff.attributes["feature"] = "SO:0000853"
        gff.attributes["score"] = "."
        gff.attributes["phase"] = "."
        gff.attributes["attributes"] = {}
        gff.attributes["attributes"]["ID"] = f"{other.species}:{other.region}"
        gff.attributes["attributes"]["species"] = other.species
        gff.attributes["attributes"]["region"] = other.region
        gff.attributes["attributes"]["strand"] = other.strand
        if self.parent:
            gff.attributes["attributes"]["parent"] = self.parent
        if seq:
            gff.attributes["attributes"]["seq"] = other.seq
            gff.attributes["attributes"]["base_seq"] = base.seq
        return gff

    @property
    def id(self):
        return f"{self.other.species}:{self.other.region}"

    @property
    def sam(self):
        flag = 64 & 255
        base = self.base
        other = self.other
        if other.strand == "-":
            flag = flag | 16
        mapq = 255
        cigar = self.cigar
        seq = other.seq.replace("-", "")
        tlen = len(seq)
        tags = [
            f"RG:Z:{other.species}",
            f"ch:Z:{other.chrom}",
            f"st:i:{other.start}",
            f"en:i:{other.end}",
            f"pt:Z:{self.parent}",
        ]
        tagstr = "\t".join(tags)
        # ret = f"{self.id}:{self.parent}\t{flag}\t{base.chrom}\t{base.start}\t{mapq}\t{cigar}\t*\t*\t{tlen}\t{seq}\t*\t{tagstr}"
        ret = f"{self.id}:{self.parent}\t{flag}\t{base.chrom}\t{base.start}\t{mapq}\t{cigar}\t*\t0\t{tlen}\t{seq}\t*\t{tagstr}"
        return ret

    @property
    def cigar(self):
        base = self.base.seq
        other = self.other.seq
        current = ""
        last = ""
        current_len = 0
        cg = ""
        for i in range(0, len(base)):
            b = base[i]
            o = other[i]
            if b == "-" and o == "-":
                raise f"Unexpected double gap in  {self.id} {i}:\n{base}\n{other}"
            elif b == o:
                current = "M"
            elif b == "-":
                current = "I"
            elif o == "-":
                current = "D"
            elif b != o:
                current = "M"

            if current != last and current_len > 0:
                cg = f"{cg}{current_len}{last}"
                current_len = 0
            current_len += 1
            last = current
        if current_len > 0:
            cg = f"{cg}{current_len}{last}"
        return cg


class EnsemblPairwiseAlignments:
    def __init__(self, response, base="triticum_aestivum"):
        self.response = response
        self.base = base
        self.__alns = list(
            map(
                lambda aln: EnsemblPairwiseAlignment(aln, base=self.base), self.response
            )
        )

    @property
    def longest(self):
        ret = self.__alns[0]
        # print("Finding longest...")
        for aln in self.alns:
            # pprint.pp(aln)
            # print(aln.gff)
            # print(f"Comparing {len(aln)} > {len(ret)}")
            if len(aln) > len(ret):
                # print("longer....")
                ret = aln
        # print(f"Returning... {ret}")
        return ret

    @property
    def alns(self):
        return self.__alns

    @property
    def gff(self):
        return map(lambda aln: aln.gff, self.alns)
