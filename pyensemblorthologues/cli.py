"""Console script for pyensemblorthologues."""

import fire


def help():
    print("pyensemblorthologues")
    print("=" * len("pyensemblorthologues"))
    print("Tool to download ortologue genes from ensembl compara")


def main():
    fire.Fire({"help": help})


if __name__ == "__main__":
    main()  # pragma: no cover
