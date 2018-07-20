#! /usr/bin/env python

from cyvcf2 import VCF
from collections import defaultdict
from argparse import ArgumentParser
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    args = get_args()
    len_dict = defaultdict(list)
    for v in VCF(args.vcf):
        if not v.INFO.get('SVTYPE') == 'TRA':
            len_dict[v.INFO.get('SVTYPE')].append(abs(v.INFO.get('SVLEN')))
    make_plot(dict_of_lengths=len_dict,
              output=args.output)


def make_plot(dict_of_lengths, output):
    standard_order = ['DEL', 'INS', 'DUP', 'INV']
    spec_order = sorted([i for i in dict_of_lengths.keys() if i not in standard_order])
    sorter = standard_order + spec_order
    names, lengths = zip(
        *sorted([(k, v) for k, v in dict_of_lengths.items()], key=lambda x: sorter.index(x[0])))
    plt.subplot(2, 1, 1)
    plt.hist(x=lengths,
             bins=[i for i in range(0, 2000, 10)],
             stacked=True,
             histtype='bar',
             label=names)
    plt.xlabel('Lenghth of structural variant')
    plt.ylabel('Number of variants')
    plt.legend(frameon=False,
               fontsize="small")

    plt.subplot(2, 1, 2)
    plt.hist(x=lengths,
             bins=[i for i in range(0, 20000, 100)],
             stacked=True,
             histtype='bar',
             label=names,
             log=True)
    plt.xlabel('Lenghth of structural variant')
    plt.ylabel('Number of variants')
    plt.legend(frameon=False,
               fontsize="small")
    plt.tight_layout()
    plt.savefig(output)


def get_args():
    parser = ArgumentParser(description="create stacked bar plot of the SV lengths split by type")
    parser.add_argument("vcf", help="vcf file to parse")
    parser.add_argument("-o", "--output", help="output file to write to", default="SV-length.png")
    return parser.parse_args()


if __name__ == '__main__':
    main()
