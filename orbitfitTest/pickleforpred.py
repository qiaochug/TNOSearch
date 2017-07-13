#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thur July 13 09:40:35 2017

@author:qiaochug
"""

import sys
import numpy as np
import argparse
import pandas as pd
import pickle
from LinkerLib import Detection
from LinkerLib import fakeObj
from math import degrees

    
def trimObj(fakelist, outfile):
    #organize the detections within each fakeObj
    #divide into four campaigns y1 = (56505,56717), y2 = (56870,57082)
    #y3 = (57235,57448), y4 = (57601,57813)

    for obj in fakelist:
        #first initiate all four campaigns
        campaign1 = []
        campaign2 = []
        campaign3 = []
        campaign4 = []
        #put each detection into corresponding campaign
        for det in obj.listobj:
            if det.mjd > 56504 and det.mjd < 56718:
                campaign1.append(det)
            elif det.mjd > 56869 and det.mjd < 57083:
                campaign2.append(det)
            elif det.mjd > 57234 and det.mjd < 57449:
                campaign3.append(det)
            elif det.mjd > 57600 and det.mjd < 57814:
                campaign4.append(det)

        #if less than 4, do not add this campaign to the campaigns list of this obj
        for campaign in (campaign1,campaign2,campaign3,campaign4):
            if len(campaign) > 3:
                obj.campaigns.append(campaign)

        #if the obj does not have any campaign with 4 or more detections
        #remove the obj from the fakelist
        if len(obj.campaigns) == 0:
            fakelist.remove(obj)
            
    pickle.dump(fakelist, (open(outfile,"wb")))
    
    return fakelist
        
def parseFile(csvName):
    df = pd.read_csv(csvName)
    fakeList = []

    #for percentage display
    listLen = len(df.index)
    counter = 0
    
    for index,row in df.iterrows():
        det = Detection(objid=row['objid'], fakeid=row['fakeid'], ra=degrees(row['ra']), dec=degrees(row['dec']), mjd=row['mjd_obs'],\
                        flux=0, expnum=row['expnum'], ccd=row['ccdnum'], band=row['band'],lookAhead=0)
        if len(fakeList) == 0 or fakeList[-1].fakeid != row['fakeid']:
            newObj = fakeObj(row['fakeid'])
            newObj.listobj.append(det)
            fakeList.append(newObj)
        else:
            fakeList[-1].listobj.append(det)
        if len(fakeList) >= 100:
            break

        #print percentage
        sys.stdout.write("\rParsing the File %2f%%"%(float(counter)/listLen*100))
        sys.stdout.flush()

    return fakeList      
    
def main():
    """
    takes in a csv file with fakes info
    the file is sorted to list the detections of the same fakeid together
    mjd from early to late
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfile', nargs = 1, default = None,
                         help = 'path to csv file')
    parser.add_argument('outputfile', nargs = 1,
                        default = None,
                        help = 'outputFile for writing the triplets')
    
    args = parser.parse_args()
    fakelist = parseFile(args.inputfile[0])
    processedfakelist = trimObj(fakelist,args.outputfile[0])

if __name__ == '__main__':
    main()
