"""Console script for pyensemblorthologues."""
import Mikado.loci
import Mikado.parsers
from Bio import SeqIO
from fire import Fire
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


class MainProgram:
    def __init__(self) -> None:

        pass


def main():
    # Fire({"help": help})

    method = "LASTZ_NET"
    species = "triticum_aestivum"
    # target_species = "triticum_turgidum"
    interval = "3B:684798558-684799943"
    compara = "plants"
    server = "http://rest.ensembl.org"
    cc = ComparaConsumer(server=server, compara=compara)
    ss = cc.species_sets(method=method, species=species)
    print(ss)
    ort = cc.regions(method=method, species=species, interval=interval)
    print(ort)

    msa = MSARegion(ort)
    print(msa.aligned())
    SeqIO.write(msa.aligned(), "example.fasta", "fasta")


if __name__ == "__main__":
    main()  # pragma: no cover
