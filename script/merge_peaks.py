#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import re
import sys
import time
import argparse
from multiprocessing import Pool
from bx.intervals.intersection import IntervalTree
from functools import cmp_to_key


parser = argparse.ArgumentParser(
    description="This script is used to get dupulicate peaks in some samples.",
    formatter_class=argparse.RawTextHelpFormatter,
    usage=(
        "python merge_peaks.py -AP peakfiles.txt -CL hg19.fa.auto.fai -O Dup_peaks.xls"
    ),
)
parser.add_argument(
    "-AP",
    "--allpeaks",
    dest="allpeaks",
    type=str,
    default=None,
    help="The all peaks files path!",
    required=True,
)
parser.add_argument(
    "-CL",
    "--chrlist",
    dest="chrlist",
    type=str,
    default=None,
    help="The useful fai file path!",
    required=True,
)
parser.add_argument(
    "-T",
    "--threads",
    dest="threads",
    type=int,
    default=24,
    help="The number of threads",
    required=False,
)
parser.add_argument(
    "-O",
    "--output",
    dest="output",
    type=str,
    default=None,
    help="The result peak file path!",
    required=True,
)


def parse_peak(file):
    files = []
    samples = []
    with open(file, "r") as f:
        for line in f:
            files.append(line.strip())
            samples.append(line.strip().split("/")[-1].split("_peaks")[0])
    return files, samples


def chromn_cmp(x: str, y: str) -> int:
    xs = re.findall("[A-z]", x)
    ys = re.findall("[A-z]", y)
    xl = len(xs)
    yl = len(ys)
    if (xl == 0) and (yl > 0):
        return -1
    if (xl > 0) and (yl == 0):
        return 1
    if (xl == 0) and (yl == 0):
        x = int(x)
        y = int(y)
        return 1 if x > y else -1 if x < y else 0
    xlp = x.rfind(xs[-1])
    ylp = y.rfind(ys[-1])
    oxl = len(x)
    oyl = len(y)
    if (xlp == (oxl - 1)) and (ylp == (oyl - 1)):
        x = x.lower()
        y = y.lower()
        return 1 if x > y else -1 if x < y else 0
    if (xlp < (oxl - 1)) and (ylp == (oyl - 1)):
        return -1
    if (xlp == (oxl - 1)) and (ylp < (oyl - 1)):
        return 1
    strx = x[:xlp].lower()
    numx = float(x[xlp + 1:])
    stry = y[:ylp].lower()
    numy = float(y[ylp + 1:])
    if strx == stry:
        return 1 if numx > numy else -1 if numx < numy else 0
    return 1 if strx > stry else -1 if strx < stry else 0


def parse_chr(file):
    chroms = []
    with open(file, "r") as f:
        for line in f:
            chroms.append(line.strip().split("\t")[0])
    chroms.sort(key=cmp_to_key(chromn_cmp))
    return chroms


def get_peak_info(peakfiles, samples, chrid):
    pkinfo = dict()
    for i in range(len(peakfiles)):
        sample = samples[i]
        f = open(peakfiles[i])
        next(f)
        for line in f:
            lines = line.strip().split("\t")
            chrom = lines[0]
            if chrom == chrid:
                peakid = lines[3]
                start = int(lines[1])
                end = int(lines[2])
                pileup = lines[-1]
                pvalue = float(lines[-3])
                infos = [chrom, start, end, pileup, pvalue, sample, peakid]
                pkinfo[peakid] = infos
        f.close()
    return pkinfo


def search_best_peak(pkinfo, tree):
    local_best_peaks = set()
    all_peak_res = dict()
    for peakid, info in pkinfo.items():
        peak_res = []
        start = info[1]
        end = info[2]
        all_res = tree.find(start, end)
        for res in all_res:
            if abs(res[1] + res[2] - start - end) < 1000:
                peak_res.append(res)
        peak_res.sort(key=lambda x: x[4])
        all_peak_res[peakid] = peak_res
        local_best_peak = peak_res[0]
        if peakid == local_best_peak[-1]:
            local_best_peaks.add(peakid)
    return local_best_peaks, all_peak_res


def search_dup_peaks(peakfiles, samples, chrid):
    out_info = []
    tree = IntervalTree()
    pkinfo = get_peak_info(peakfiles, samples, chrid)
    for peakid, info in pkinfo.items():
        tree.insert(info[1], info[2], info)
    local_best_peaks, all_peak_res = search_best_peak(pkinfo, tree)
    del pkinfo
    del tree
    for pkid in local_best_peaks:
        from_samples = set()
        for i in all_peak_res[pkid]:
            from_samples.add(i[-2])
        duplicates = len(from_samples)
        if duplicates < 2:
            continue
        pk_start = all_peak_res[pkid][0][1]
        pk_end = all_peak_res[pkid][0][2]
        from_samples = ",".join(list(from_samples))
        row_info = [
            pkid,
            chrid,
            str(pk_start),
            str(pk_end),
            str(duplicates),
            from_samples
        ]
        for sample in samples:
            if sample not in from_samples:
                df_pileup = "0"
            else:
                for peak in all_peak_res[pkid]:
                    if sample == peak[-2]:
                        df_pileup = peak[3]
            row_info.append(df_pileup)
        out_info.append(row_info)
    out_info.sort(key=lambda x: (int(x[2]), int(x[3])))
    return out_info


def parallel_process(peakfiles, samples, chroms, output, threads):
    column_names = [
        "PeakID",
        "Peak_chr",
        "Peak_Start",
        "Peak_end",
        "Dup_num",
        "Dup_samples",
    ]
    column_names += samples
    out = open(output, "w")
    out.write("\t".join(column_names) + "\n")
    workers_thread = []
    pool = Pool(processes=threads)
    for chrid in chroms:
        w = pool.apply_async(
            search_dup_peaks,
            (peakfiles, samples, chrid),
        )
        workers_thread.append(w)
    pool.close()
    pool.join()
    for w in workers_thread:
        result = w.get()
        for dup_info in result:
            out.write("\t".join(dup_info) + "\n")


def main():
    argv = vars(parser.parse_args())
    allpeaks = argv["allpeaks"]
    chrlist = argv["chrlist"]
    threads = argv["threads"]
    output = argv["output"]
    peakfiles, samples = parse_peak(allpeaks)
    chroms = parse_chr(chrlist)
    parallel_process(peakfiles, samples, chroms, output, threads)


if __name__ == "__main__":
    print(time.strftime("Start %Y-%m-%d %H:%M:%S", time.localtime()))
    main()
    print(time.strftime("Done %Y-%m-%d %H:%M:%S", time.localtime()))
