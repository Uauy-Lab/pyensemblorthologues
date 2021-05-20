import requests

"""Main module."""


class ComparaConsumer:
    def __init__(self, server="http://rest.ensembl.org", compara="plants"):
        self.server = server
        self.compara = compara

    def request(self, url, args):
        full_url = f"{self.server}/{url}?{args}"
        print(full_url)
        r = requests.get(full_url, headers={"Content-Type": "application/json"})
        if not r.ok:
            print(f"Error reading: `{full_url}`")
            # r.raise_for_status()
        # print(r.text)
        return r.json()

    def region(
        self,
        method="LASTZ_NET",
        species="triticum_aestivum",
        target_species="triticum_turgidum",
        interval="3B:684798558-684799943",
    ):
        url = f"alignment/region/{species}/{interval}"
        args = f"compara={self.compara};method={method};species_set={species};species_set={target_species}"
        return self.request(url, args)

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
        return ret

        # http://rest.ensembl.org/info/compara/species_sets/LASTZ_NET?content-type=application/json;compara=plants;

    def regions(
        self,
        mathod="LASTZ_NET",
        species="triticum_aestivum",
        interval="3B:684798558-684799943",
    ):
        species_sets = self.species_sets(method=mathod, species=species)
        for sp in species_sets:
            alignment = self.region(
                species=species, target_species=sp, interval=interval
            )
            print(alignment)


if __name__ == "__main__":
    cc = ComparaConsumer()
    ss = cc.species_sets()
    ort = cc.regions()

    print(ss)
    print(ort)
