"""Console script for pyensemblorthologues."""
import pprint

import Mikado.loci
import Mikado.parsers
from Bio import SeqIO
from fire import Fire
from Mikado.parsers.GFF import GffLine
from pyensemblorthologues.MSARegions import MSARegion

from .compara_consumer import ComparaConsumer


def strip_utr(gene: Mikado.loci.Gene):

    for transcript in gene:
        assert isinstance(transcript, Mikado.loci.Transcript)
        if len(transcript.combined_cds) > 0:
            transcript.exons = transcript.combined_cds.copy()
            transcript.start = min([_[0] for _ in transcript.combined_cds])
            transcript.end = max([_[1] for _ in transcript.combined_cds])
            transcript.combined_utr = []
        transcript.finalize()
        # transcript.remove_utrs()

    gene.finalize()
    return gene


def help():
    print("pyensemblorthologues")
    print("=" * len("pyensemblorthologues"))
    print("Tool to download ortologue genes from ensembl compara")


def region_for_gene(row, flank=5000):
    id = row.id.replace("gene:", "")
    row.id = f"{id}:{flank}"
    start = row.start - flank
    end = row.end + flank
    row.feature = "SO:0000239"
    row.start = start
    row.end = end

    return f"{row.chrom}:{start}-{end}"


class Compara:
    def extract(
        self,
        gff,
        method="LASTZ_NET",
        species="triticum_aestivum",
        flank=2000,
        compara="plants",
        server="http://rest.ensembl.org",
        output="output.gff",
        longest=False,
        sequence=False,
    ):
        cc = ComparaConsumer(server=server, compara=compara)
        # print(gff)
        # Examples to try: TraesCS7A02G175200 (Nikolai) TraesCS6A02G313800 (Andy)
        parser = Mikado.parsers.parser_factory(gff, "gff3")
        i = 0
        f = open(output, "w")
        for row in parser:
            if row.is_gene is True and row.attributes["biotype"] == "protein_coding":
                interval = region_for_gene(row, flank=flank)
                print(interval)
                print(row)
                print(row.attributes)
                ort = cc.regions(
                    method=method, species=species, interval=interval, longest=longest, parent=row.id
                )
                f.write(str(row))
                f.write("\n")
                for aln in ort:
                    f.write(GffLine.string_from_dict(aln.gff(seq=sequence).attributes))
                    f.write("\n")
                i += 1
            if i > 10:
                break
        f.close()
        # ss = cc.species_sets(method=method, species=species)

    def msa(self, ort):
        if len(ort) < 5:
            pass  # should be a continue.
        msa = MSARegion(ort)
        print(msa.unaligned)
        # print(msa.aligned())
        SeqIO.write(msa.aligned(), f"{id}.fasta", "fasta")
        pass


class Pipeline:
    def __init__(self) -> None:
        self.compara = Compara()


def main():
    # Fire({"help": help})
    Fire(Pipeline)


def old_main():
    method = "LASTZ_NET"
    species = "triticum_aestivum"
    # target_species = "triticum_turgidum"
    # interval = "3B:684798558-684799943"
    flank = 2000
    start = 575672636
    end = 575673382
    start = start - flank
    end = start - flank
    chr = "1A"
    interval = f"{chr}:{start}-{end}"
    compara = "plants"
    server = "http://rest.ensembl.org"
    cc = ComparaConsumer(server=server, compara=compara)
    # ss = cc.species_sets(method=method, species=species)

    ort = cc.regions(method=method, species=species, interval=interval)
    print(ort)

    msa = MSARegion(ort)
    print(msa.aligned())
    SeqIO.write(msa.aligned(), "example.fasta", "fasta")


if __name__ == "__main__":
    main()  # pragma: no cover
