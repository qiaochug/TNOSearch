#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This program add the MJD feature to the training set
The MJD feature together with the dvlon dvlat cos12 cos23
should be quiet sufficient for a human expert
to tell whether the triplet is good

Created on Thu Jun 22 14:50:22 2017

@author: qiaochug
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def addmjd(trainingset, lookuptable):
    ts = pd.read_csv(trainingset)
    lookup = pd.read_csv(lookuptable)
    ts = ts.merge(lookup, left_on='id1',right_on='SNOBJID',how = 'left',sort = False)
    ts = ts.merge(lookup, left_on='id2',right_on='SNOBJID',how = 'left',sort = False)
    ts = ts.merge(lookup, left_on='id3',right_on='SNOBJID',how = 'left',sort = False)
        
    return ts

def newtrainingset(ts, output):
    ts.to_csv(output)
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('trainingset', nargs = 1, default = None,
                         help = 'path to trainingset without mjd')
    parser.add_argument('lookuptable', nargs = 1, default = None,
                         help = 'path to table containing the mjd lookup')
    parser.add_argument('completetrainingset', nargs = 1, default = None,
                         help = 'path to output the completed training set with mjd')
    
    args = parser.parse_args()
    
    ts = addmjd(args.trainingset[0], args.lookuptable[0])
    newtrainingset(ts, args.completetrainingset[0])
    
if __name__ == '__main__':
        main()