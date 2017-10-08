#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 15:27:50 2017
vcf2able.py
@author: scott
"""

import argparse
from collections import defaultdict
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-i', "--vcfin", type=str, required=True,
                    help="vcffile")
parser.add_argument('-b', "--block_size", type=int, required=True,
                    help="size of linkage blocks")
parser.add_argument('-d', "--pedfile", type=str, required=True,
                    help="path to pedfile with pop and individuals")
parser.add_argument('-p', "--pop", type=str, nargs='+', required=False,
                    default='', help="which pops, if blank uses all")
parser.add_argument('-s', "--size", type=str, nargs='+', required=False,
                    default='', help="how many from each pop, blank uses all")
parser.add_argument("--bin", action="store_true",
                    help="create pseudo-ms with binary instead of gt")
args = parser.parse_args()

# TODO: make binary parser for bin option


def firstchrom(vcfin):
    """
    """
    with open(vcfin, 'r') as vcf:
        for line in vcf:
            if line.startswith("#CHROM"):
                pop_iix = line.strip().split()
                line = vcf.next()
                chrom = line.strip().split()[0]
                pos1 = line.strip().split()[1]
                break
    return(chrom, int(pos1), pop_iix)


def makeable(i, pop, p, x, abledict):
    """
    """
    ra = x[3]
    aa = x[4]
    gt = x[p].split(":")[0]
    if gt.count("0") == 2:
        abledict[pop][2*i] += ra
        abledict[pop][2*i+1] += ra
    elif gt.count("1") == 2:
        abledict[pop][2*i] += aa
        abledict[pop][2*i+1] += aa
    else:
        if gt == "0/1" or gt == "0|1":
            abledict[pop][2*i] += ra
            abledict[pop][2*i+1] += aa
        elif gt == "1/0" or gt == "1|0":
            abledict[pop][2*i] += aa
            abledict[pop][2*i+1] += ra
        else:
            abledict[pop][2*i] += "N"
            abledict[pop][2*i+1] += "N"
    return(abledict)


def able2vcf(vcfin, block, popinfo, pops, sizes, rand=True):
    """
    """
    # parse ped
    peddict = defaultdict(list)
    with open(popinfo, 'r') as ped:
        for line in ped:
            if line.strip():
                x = line.strip().split()
                if (x[0] in pops) or (pops == '0'):
                    peddict[x[0]].append(x[1])
            else:
                continue
    poplist = pops
    if sizes:
        sizes = map(int, sizes)
        if rand:
            for pop in poplist:
                i = pops.index(pop)
                peddict[pop] = np.random.choice(peddict[pop], sizes[i],
                                                replace=False)
        else:
            for pop in poplist:
                i = pops.index(pop)
                peddict[pop] = peddict[pop][:sizes[i]]
    # get initial
    chrom, pos1, pop_iix = firstchrom(vcfin)
    # make able file
    abledict = {}
    for pop in poplist:
        abledict[pop] = [''] * 2 * len(peddict[pop])
    b = 0
    with open("able.in", 'w') as able:
        with open(vcfin, 'r') as vcf:
            for line in vcf:
                if not line.startswith("#"):
                    x = line.strip().split()
                    pos = int(x[1])
                    if (pos <= (block + b)) and (x[0] == chrom):
                        for pop in poplist:
                            for i, sample in enumerate(peddict[pop]):
                                p = pop_iix.index(sample)
                                abledict = makeable(i, pop, p, x, abledict)
                    else:
                        able.write("\n//\n#{}_{}-{}\n".format(chrom, b,
                                                              block+b))
                        if abledict[pop][0] == '':
                            pass
                        else:
                            for pop in poplist:
                                for samp in abledict[pop]:
                                    able.write("{}\n".format(samp))
                        b += block
                        for pop in poplist:
                            abledict[pop] = [''] * 2 * len(peddict[pop])
                        # restart
                        if (pos <= (block + b)) and (x[0] == chrom):
                            for pop in poplist:
                                abledict[pop] = [''] * 2 * len(peddict[pop])
                                for i, sample in enumerate(peddict[pop]):
                                    p = pop_iix.index(sample)
                                    abledict = makeable(i, pop, p, x, abledict)
                        else:
                            able.write("\n//\n#{}_{}-{}\n".format(chrom, b,
                                                                  block+b))
                            b += block
                    # update chrom
                    chrom = x[0]


def able2vcf_bin(vcf, block, popinfo, pops):
    """
    """


if __name__ == "__main__":
    vcf = args.vcfin
    block = int(args.block_size)
    popinfo = args.pedfile
    pops = args.pop
    sizes = args.size
    if args.bin:
        able2vcf_bin(vcf, block, popinfo, pops, sizes)
    else:
        able2vcf(vcf, block, popinfo, pops, sizes)
