#!/usr/bin/env python

from smact.properties import band_gap_Harrison

# python band_gap_simple.py --h
   # usage: band_gap_simple.py [-h] [-a ANION] [-c CATION] [-d DISTANCE] [-v]

#python band_gap_simple.py -c Mg -a Cl -d 2.38
   # 3.8944137939094166

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Compound band gap estimates from elemental data."
    )
    parser.add_argument(
        "-a", "--anion", type=str, help="Element symbol for anion."
    )
    parser.add_argument(
        "-c", "--cation", type=str, help="Element symbol for cation."
    )
    parser.add_argument(
        "-d", "--distance", type=float, help="Internuclear separation."
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="More Verbose output."
    )
    args = parser.parse_args()

    print(
        band_gap_Harrison(
            verbose=args.verbose,
            anion=args.anion,
            cation=args.cation,
            distance=args.distance,
        )
    )
