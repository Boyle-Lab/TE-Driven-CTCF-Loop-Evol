#!/usr/bin/python

import sys, re, os
from optparse import OptionParser
import pandas as pd
import numpy as np
import tabix as tb


"""
Subroutine definitions
"""

def trim_anchors(row, peaks, opt):
    """
    Given a row from the 3D contacts dataframe, trim both anchors to the
    boundaries of the (best) overlapping ChIP-seq peak.
    """
    left_trimmed = right_trimmed = False
    
    peak1 = find_best_peak(row[0],
                           row[1],
                           row[2],
                           "proximal",
                           peaks,
                           opt)
    peak2 = find_best_peak(row[3],
                           row[4],
                           row[5],
                           "distal",
                           peaks,
                           opt)
    if peak1.shape[0] > 0:
        # Trim anchors and add summit1
        row['start1'] = int(peak1.chromStart)
        row['end1'] = int(peak1.chromEnd)
        row['summit1'] = int(peak1.peak) + int(peak1.chromStart)
        left_trimmed = True
    else:
        row['summit1'] = int( (row['start1'] + row['end2']) / 2)
    if peak2.shape[0] > 0:
        # Trim anchors and add summit2
        row['start2'] = int(peak2.chromStart)
        row['end2'] = int(peak2.chromEnd)
        row['summit2'] = int(peak2.peak) + int(peak2.chromStart)
        right_trimmed = True
    else:
        row['summit2'] = int( (row['start2'] + row['end2']) / 2)
    return row, left_trimmed, right_trimmed


def find_best_peak(chrom, chromStart, chromEnd, end, peaks, opt):
    """
    Given a genomic interval, see if any ChIP-seq peaks overlap and
    return the row corresponding to the best overlapping peak, or
    an empty dataframe if none overlap.
    """
    feats = get_features(chrom,
                         chromStart,
                         chromEnd,
                         peaks,
                         ["chrom",
                          "chromStart",
                          "chromEnd",
                          "name",
                          "score",
                          "strand",
                          "signalValue",
                          "pValue",
                          "qValue",
                          "peak"],
                         ["string",
                          "int64",
                          "int64",
                          "string",
                          "int64",
                          "string",
                          "float64",
                          "float64",
                          "float64",
                          "int64"])
    return choose_feat(feats, "signalValue", end)


def choose_feat(feats, col, end):
    """
    Given a set of features, choose the one with maximum value.
    In case of ties, use the closest to the proximal or distal
    end, depending on whether this is an upstream or downstream
    anchor.
    """
    if feats.shape[0] == 0:
        return feats
    f = feats[feats[col] == feats[col].max()]
    if f.shape[0] > 1:
        if end == "proximal":
            return f[f.peak == f.peak.max()]
        else:
            return f[f.peak == f.peak.min()]
    else:
        return f


def get_features(chrom, start, end, feats_f, names, dtypes):
    """
    Get a pandas dataframe of features within a given anchor. Uses Tabix.
    """
    if len(names) != len(dtypes):
        sys.stderr.write("get_features: ERROR -- names and dtypes must have same length!\n")
        exit(1)
    feats = feats_f.querys( '{}:{}-{}'.format(chrom,
                                              start,
                                              end) )
    # Convert to a pandas dataframe with the supplied column names
    feats = pd.DataFrame(list(feats), columns = names)
    # Convert datatypes for any numeric columns (strings should be fine as is)
    for i in range(0, len(names)):
        d = dtypes[i]
        if d == 'int' or d == 'int64' or d == 'float' or d == 'float64':
            feats[names[i]] = pd.to_numeric(feats[names[i]])
    return feats
        
"""
Main Definitions
"""

if __name__ == "__main__":
    parser = OptionParser(description='Given a set of 3D contacts from a ChIA-PET experiment, in bedpe format, find the highest-scoring ChIP-seq peaks in each anchor and trim the anchors to correspond to the boundaries of this element. Store the summit of the element as an additional column (summit1 and summit2).')
    parser.add_option('-b', '--bedpe_file', action="store", type="string", dest="bedpe", metavar="<file>",
                        help='Path to the 3D contacts file, in bedpe format. E.g., from MANGO predictions.')
    parser.add_option('-c', '--ctcf_chipseq_file', action="store", type="string", dest="peaks", metavar="<file>",
                        help='Path to CTCF ChIP-seq data file. See help document for format requirements.')
    parser.add_option('-l', '--min_loop_size', type="int", default=0,
                        help='Minimum distance between CTCF summits to retain a loop.')
    parser.add_option('-u', '--max_loop_size', type="int", default=0,
                        help='Maximum distance between CTCF summits to retain a loop.')
    parser.add_option('-k', '--keep_unoccupied', action="store_true", default=False,
                        help='Keep loops in which one or both anchors do not contain a CTCF peak.')

    (opt, args) = parser.parse_args(sys.argv)


    """
    Load the bedpe file into a pandas dataframe
    """
    loops = pd.read_table(opt.bedpe)
    
    """
    Go through the bedpe table row by row, checking each loop
    for ChIP-seq overlaps. If overlap(s) are found, trim the
    anchor(s) to the boundaries of the (best) peak in each anchor.
    """
    peaks = tb.open(opt.peaks)
    for idx, row in loops.iterrows():
        new_row, left_trimmed, right_trimmed = trim_anchors(row, peaks, opt)
        good = False
        if opt.keep_unoccupied:
            good = True
        else:
            if left_trimmed and right_trimmed:
                good = True
        if good:
            sys.stdout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(new_row.chrom1,
                                                                       new_row.start1,
                                                                       new_row.end1,
                                                                       new_row.chrom2,
                                                                       new_row.start2,
                                                                       new_row.end2,
                                                                       new_row.IAB,
                                                                       new_row.FDR))
