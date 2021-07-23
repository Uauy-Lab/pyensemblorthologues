import os
import re
from os.path import exists

import psutil
import pyensemblorthologues
from pyfaidx import Fasta


def get_pt_from_sam_line(line):

    tokens = line.rstrip().split("\t")
    ret = None
    for tok in tokens:
        if tok.startswith("pt:Z:"):
            ret = re.sub("pt:Z:", "", tok)
    return ret


def cmdLine():
    my_process = psutil.Process(os.getpid())
    return my_process.cmdline()


class ComparaWriter:
    def __init__(
        self, path, compara_consumer=None, format="sam", reference=None
    ) -> None:
        self.path = path
        self.cc = compara_consumer
        self.reference = reference
        self.format = format
        self.file = None
        self.last_record_id = self.__last_record_id()

    def __last_record_id(self):
        if not exists(self.path):
            return None
        with open(self.path, "rb") as f:
            f.seek(-2, os.SEEK_END)
            while f.read(1) != b"\n":
                f.seek(-2, os.SEEK_CUR)
            last_line = f.readline().decode()
            if last_line.startswith("@"):
                return None
        return get_pt_from_sam_line(last_line)

    def header(self, species, method):
        header = list()
        if self.reference is not None:
            ref = Fasta(self.reference)
            recs = list(
                map(lambda record: f"@SQ\tSN:{record.name}\tLN:{len(record)}", ref)
            )
            header.extend(recs)
        header.extend(self.cc.sam_header(method=method, species=species))

        cmd = " ".join(cmdLine())
        # print(cmd)
        header.append(
            f"@PG\tID:pyensemblorthologues\tPN:pyensemblorthologues\tVN:{pyensemblorthologues.__version__}\tCL:{cmd}"
        )
        return header

    def close(self):
        if self.file is not None:
            self.file.close()
        self.file = None

    def open(self, species="triticum_aestivum", method="LASTZ_NET"):
        if self.file is not None:
            return
        if self.last_record_id is None:
            self.file = open(self.path, "w")
            header = self.header(species, method)
            # print(header)
            self.file.write("\n".join(header))
            self.file.write("\n")
            self.file.close()

        self.file = open(self.path, "a")

    def write_aln(self, aln):
        sam_line = aln.sam
        self.last_record_id = get_pt_from_sam_line(sam_line)
        self.file.write(sam_line)
        self.file.write("\n")
