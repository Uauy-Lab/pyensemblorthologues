import json
import time
from os import wait
from pprint import pprint
from typing import Any

import requests
from pyensemblorthologues.ensembl_pairwise_alignment import EnsemblPairwiseAlignments


class ComparaConsumer:
    def __init__(self, server="http://rest.ensembl.org", compara="plants"):
        self.server = server
        self.compara = compara

    def request(self, url, args, max_try=10, max_wait=600) -> Any or None:
        full_url = f"{self.server}/{url}?{args}"
        # print(full_url)
        r = requests.get(full_url, headers={"Content-Type": "application/json"})
        wait_seconds = 1
        current_try = 0
        while current_try < max_try:
            if r.ok:
                return r.json()
            if r.status_code == 400:
                return None
            current_try = current_try + 1
            print(
                f"Waiting {wait_seconds} seconds to retry({current_try})  #{full_url}"
            )
            time.sleep(wait_seconds)
            wait_seconds = wait_seconds * 2
            if wait_seconds > max_wait:
                wait_seconds = max_wait
        r.raise_for_status()

    def region(
        self,
        method="LASTZ_NET",
        species="triticum_aestivum",
        target_species="triticum_turgidum",
        interval="3B:684798558-684799943",
    ):

        url = f"alignment/region/{species}/{interval}"
        args = f"compara={self.compara};method={method};species_set={species};species_set={target_species}"
        ret = None
        try
            ret = self.request(url, args)
        except 
            with open("errors.txt","a") as out :
                out.write(f"Failed:{species}:{target}:{interval}:{method}")
        return ret 


    def species_sets(self, method="LASTZ_NET", species="triticum_aestivum"):
        url = f"info/compara/species_sets/{method}"
        args = f"compara={self.compara}"
        ss = self.request(url, args)
        ret = []
        for group in ss:
            species_set = group["species_set"]
            sp = list(filter(lambda s: s != species, species_set))
            if len(sp) == 1:
                ret.append(sp[0])
        return list(set(ret))

    def regions(
        self,
        method="LASTZ_NET",
        species="triticum_aestivum",
        interval="3B:684798558-684799943",
        longest=False,
        parent=None,
    ):
        species_sets = self.species_sets(method=method, species=species)
        ret = []
        for sp in species_sets:
            alignment = self.region(
                species=species, target_species=sp, interval=interval
            )
            if not alignment:
                continue
            epas = EnsemblPairwiseAlignments(alignment, base=species)
            if longest:
                ret.append(epas.longest)
            else:
                ret.extend(epas.alns)
        for aln in ret:
            aln.parent = parent
        return ret

    def sam_header(
        self, method="LASTZ_NET", species="triticum_aestivum", reference=None
    ):
        species_sets = self.species_sets(method=method, species=species)
        lines = []

        for species in species_sets:
            lines.append(f"@RG\tID:{species}\tSM:{species}")

        return lines
