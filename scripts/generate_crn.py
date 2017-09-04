#!/usr/bin/env python

import random
import argparse


def rndSpecies(n, prefix='X'):
    return prefix + str(random.choice(range(1, n + 1)))


def reactand(s, n):
    n = random.choice(range(1, n + 1))
    sp = []
    for i in range(n):
        sp.append(rndSpecies(s))
    return ' + '.join(sp)


def arrow(no_reversible):
    if no_reversible or random.choice([1, 2]) == 2:
        return '->'
    return '<=>'


def main():
    """ Generate a random CRN """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-s", "--species", type=int, default=5,
                        help="Specify a number of Species")
    parser.add_argument("-r", "--reactions", type=int, default=10,
                        help="Specify a number of reactions")
    parser.add_argument("-e", "--maxeducts", type=int, default=2,
                        help="Specify a maximum number of reaction educts")
    parser.add_argument("-p", "--maxproducts", type=int, default=2,
                        help="Specify a maximum number of reaction products")
    parser.add_argument("--no_reversible", action="store_true",
                        help="Forbid reversible reactions")
    args = parser.parse_args()

    for i in range(0, args.reactions):
        print "{} {} {}".format(
            reactand(args.species, args.maxeducts),
            arrow(args.no_reversible),
            reactand(args.species, args.maxproducts))


if __name__ == '__main__':
    main()
