"""Console script for pyensemblorthologues."""
import os
import pathlib
import pprint

import Mikado.loci
import Mikado.parsers
import psutil
from Bio import Seq, SeqIO
from fire import Fire
from Mikado.parsers.GFF import GffLine
from pyensemblorthologues.MSARegions import MSARegion
from pyfaidx import Fasta

from .compara_consumer import ComparaConsumer


def ouput_path(gff=None, output=None, chromosome=None):
    if output is None and gff is None:
        raise "Output folder for the pipline or GFF file to infer the path is required"
    if output is None:
        output = os.path.splitext(gff)[0]
    if chromosome is not None:
        output = f"{output}_{chromosome}"
    pathlib.Path(output).mkdir(parents=True, exist_ok=True)
    return output


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


def cmdLine():
    my_process = psutil.Process(os.getpid())
    return my_process.cmdline()


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
        output=None,
        longest=False,
        reference=None,
        format="sam",
        chromosome = None
    ):
        cc = ComparaConsumer(server=server, compara=compara)
        # print(gff)
        # Examples to try: TraesCS7A02G175200 (Nikolai) TraesCS6A02G313800 (Andy)
        parser = Mikado.parsers.parser_factory(gff, "gff3")
        i = 0
        output = ouput_path(gff=gff, output=output, chromosome=chromosome)
        print(f"output: {output}")

        output_aln = f"{output}/compara_aln.{format}"
        f = open(output_aln, "w")
        if format == "sam":
            header = list()

            if reference is not None:
                ref = Fasta(reference)
                recs = list(
                    map(lambda record: f"@SQ\tSN:{record.name}\tLN:{len(record)}", ref)
                )

                header.extend(recs)
            header.extend(cc.sam_header(method=method, species=species))

            cmd = " ".join(cmdLine())
            print(cmd)
            header.append(
                f"@PG\tID:pyensemblorthologues\tPN:pyensemblorthologues\tCL:{cmd}"
            )
            f.write("\n".join(header))
            f.write("\n")

        for row in parser:
            if chromosome is not None and row.chrom != chromosome:
                continue

            if row.is_gene is True and row.attributes["biotype"] == "protein_coding":
                interval = region_for_gene(row, flank=flank)
                print(interval)
                # row.attributes["species"] = species
                # print(row)
                # print(row.attributes)
                ort = cc.regions(
                    method=method,
                    species=species,
                    interval=interval,
                    longest=longest,
                    parent=row.id,
                )

                ort.sort()
                if format == "gff":
                    f.write(f"{str(row)};species={species}")
                    f.write("\n")
                for aln in ort:
                    if format == "gff":
                        f.write(GffLine.string_from_dict(aln.gff(seq=True).attributes))
                    elif format == "sam":
                        f.write(aln.sam)
                    f.write("\n")
                i += 1
            #if i > 10:
            #    break
        f.close()
        # ss = cc.species_sets(method=method, species=species)

    def msa(self, output=None, gff=None, reference=None):
        output = ouput_path(gff=gff, output=output)
        parser = Mikado.parsers.parser_factory(f"{output}/compara_aln.gff3", "gff3")
        msa = None
        current_parent = None
        if reference:
            ref = Fasta(reference)
        for row in parser:
            if row.feature == "SO:0000239":  # This is the parent region
                if msa:
                    path = f"{output}/{current_parent.id}.fasta"
                    SeqIO.write(msa.regions, path, "fasta")
                msa = MSARegion(
                    reference=reference,
                    species=row.attributes["species"],
                    regions=list(),
                )
                seq = Seq.Seq(str(ref[row.chrom][row.start - 1 : row.end]))
                msa.add_sequence(seq, row.id, row.attributes["species"])
                current_parent = row
            else:
                seq = Seq.Seq(row.attributes["seq"])
                # if row.attributes["strand" ] == "-":
                #     seq = seq.reverse_complement()
                offset = row.start - current_parent.start
                msa.add_sequence(seq, row.id, row.attributes["species"], offset=offset)
            pprint.pp(msa)
        if msa:
            path = f"{output}/{current_parent.id}.fasta"
            SeqIO.write(msa.regions, path, "fasta")
        # if len(ort) < 5:
        #     pass  # should be a continue.
        # msa = MSARegion(ort)
        # print(msa.unaligned)
        # # print(msa.aligned())
        # SeqIO.write(msa.aligned(), f"{id}.fasta", "fasta")
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
