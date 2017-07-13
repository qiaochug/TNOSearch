#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thur July 13 14:05:54 2017

@author:qiaochug
"""

import sys
import numpy as np
import argparse
import pandas as pd
import pickle
from LinkerLib import Detection
from LinkerLib import fakeObj
from LinkerLib import Triplet
from astropy.time import Time

def mjdToDate(mjd):
    dt = Time(mjd, format='mjd').datetime
    str = dt.isoformat()
    str = str.replace('-','/')
    str = str.replace('T',' ')
    return str

def compareFit(picklein, outtxt):
    fakelist = pickle.load(open(picklein,"rb"))
    total = len(fakelist)
    count = 0
    with open(outtxt, "w+") as output:
        for obj in fakelist:
            #output.write('***'+str(obj.fakeid)+'***\n')
            print(float(count)/total)
            for campaign in obj.campaigns:
                #output.write('-campaign-\n')
                trip = Triplet([campaign[0], campaign[1], campaign[2]])
                #output.write('First three:('+str(campaign[0].mjd)+','+str(campaign[0].ra)+','+str(campaign[0].dec)+'), ')
                #output.write('('+str(campaign[1].mjd)+','+str(campaign[1].ra)+','+str(campaign[1].dec)+'), ')
                #output.write('('+str(campaign[2].mjd)+','+str(campaign[2].ra)+','+str(campaign[2].dec)+')\n')
                trip.calcOrbit()
                for index in np.arange(3,len(campaign)):
                    #output.write('Det '+str(index)+': mjd '+str(campaign[index].mjd)+'\n')
                    #output.write('true: RA '+str(campaign[index].ra)+', DEC '+str(campaign[index].dec)+'\n')
                    (ra,dec),err = trip.predictPos(mjdToDate(campaign[index].mjd))
                    #output.write('pred: RA '+str(round(ra,6))+', DEC '+str(round(dec,6))+'\n')
                    output.write('                                          differ by RA '+str(round(campaign[index].ra-ra,2)) + ',DEC '+\
                                 str(round(campaign[index].dec -dec,2))+ ' in '+ str(round(campaign[index].mjd-campaign[2].mjd,2))+' days\n')
                #output.write('-----------------------------\n')
            #output.write('**********************\n\n')
            count= count+1
        

def main():
    """
    take in pickle file that contains a list of
    fake obj, which each containing a campaigns attribute with
    a list of campaigns with 4 or 5 detections within
    """
    
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfile', nargs = 1, default = None,
                         help = 'path to pickle file')
    parser.add_argument('outputfile', nargs = 1,
                        default = None,
                        help = 'outputFile for writing the prediction')
    
    args = parser.parse_args()
    compareFit(args.inputfile[0],args.outputfile[0])

if __name__ == '__main__':
    main()
