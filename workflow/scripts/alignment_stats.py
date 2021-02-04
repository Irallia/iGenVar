#! /usr/bin/env python

import pysam
from argparse import ArgumentParser


class Bamstats(object):
    def __init__(self, bam):
        samfile = pysam.AlignmentFile(bam, "rb")
        self.name = ".".join(bam.split(".")[:-1])
        self.mapped = samfile.mapped
        self.unmapped = samfile.unmapped
        self.total = samfile.mapped + samfile.unmapped
        self.mapped_p = round(self.mapped / self.total, ndigits=2)


def main():
    args = get_args()
    stats = [Bamstats(b) for b in args.bams]
    with open(args.output, 'w') as output:
        output.write("Samples:\tTotal\t{}\n".format('\t'.join([s.name for s in stats])))
        output.write("Mapped:\t{0}\t{1}\n".format(sum([s.mapped for s in stats]), '\t'.join([str(s.mapped) for s in stats])))
        output.write("Unmapped:\t{0}\t{1}\n".format(sum([s.unmapped for s in stats]), '\t'.join([str(s.unmapped) for s in stats])))
        output.write("Mapped_fraction:\t{0}\t{1}\n".format(round(sum([s.mapped for s in stats]) / sum([s.total for s in stats]), ndigits=3), '\t'.join([str(s.mapped_p) for s in stats])))


def get_args():
    parser = ArgumentParser(description="Creates basic summary statistics from a set of bam files.")
    parser.add_argument("-o", "--output",
                        default="bam_stats.txt",
                        help="path and name of output file. Directories must exist.")
    parser.add_argument("bams",
                        nargs='+',
                        help="the bam files to use for the summary")
    return parser.parse_args()


if __name__ == '__main__':
    main()
