from .compara_consumer import ComparaConsumer

"""Main module."""


if __name__ == "__main__":
    cc = ComparaConsumer()
    ss = cc.species_sets()
    ort = cc.regions()

    print(ss)
    print(ort)
