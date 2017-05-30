#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 26 13:39:09 2017

@author: qiaochug
"""

import numpy as np
import matplotlib.pyplot as plt
import conversionsNoAst as cnv
from math import radians, degrees, atan2, asin, sin, cos, pi,sqrt
from vectorCalc import predict,maxVel,maxDegPerDay
diffs = 1
def startObsDis(sNite,sRa,sDec,eNite,eRa,eDec):
    """
    Para: start and end observation
    Return: range of probably distance
    """
    global diffs
    #divided the 20-600 AU into steps of distance
    steps = 20
    ds = np.linspace(20,600,steps)
    diffs = np.zeros(steps)
    Index = 0
    for d in ds:
        #compute the vector
        (betaS, gammaS) = cnv.equatorialToEcliptic(sRa, sDec, sNite)
        gammaN, betaN = predict(np.radians(gammaS), np.radians(betaS), (eNite - sNite),d)
        (raN, decN) = cnv.eclipticToEquatorial(degrees(betaN), degrees(gammaN), eNite)
        #changable d, different predict function! note when integrating
        #print(sRa,"sDec", sDec, raN, decN)
        radius = maxDegPerDay(d)*(eNite-sNite)
        #print(radius)
        diff = sqrt((eRa-raN)**2+(eDec-decN)**2)
        #store the absolute difference for each distance in the array diffs
        # indexed same to ds
        if diff <= radius:
            diffs[Index] = 0
        else:
            diffs[Index] = abs(diff - radius)
        Index = Index + 1
    #print(diffs)
    #The diff array shows a bowl shape, truncated below 0,
    #First we find the point from which our d range should radiate from
    radP = -1
    # first step: find the smallest value in the array
    Index = 1
    minv = diffs[0]
    while Index < steps:
        if diffs[Index] < minv:
            minv = diffs[Index]
        Index = Index + 1
    #find the index of first minv walking from 0 to 19
    Index = 0
    firstZInd = -1
    while Index < steps:
        if diffs[Index] == minv:
            firstZInd = Index
            break
        Index = Index + 1
    #find the index of last minv
    Index = steps - 1
    lastZInd = -1
    while Index > -1:
        if diffs[Index] == minv:
            lastZInd = Index
            break
        Index = Index - 1
    # if both of the first and last index of minv are same there is one
    # minv in the array, if they are not equal but not -1, there are at least two
    if firstZInd == lastZInd:
        radP = firstZInd
    else:
        radP = round((firstZInd + lastZInd)/2)
    # take 5 values around the radiation points as range of d
    # consider the edge cases where the radP is smaller than 2 or larger than steps - 2
    range = ([ds[max(0,radP-2)],ds[min(radP+2, steps -1)]])      
    return range
    
def startObsPoint(sObsDis, sNite, sRa,sDec):
    """
    Para: start observation, and suspected distance(one)
    Return: A point in heliocentric coordinates
    """
    
def endObsPoint(startP, eNite,eRa,eDec):
    """
    Para: start point in heliocentric coordinates, end observation
    Return: A point in heliocentric coordinates
    """
    
def iFinding(startP, endP):
    """
    Para: start point and end point in heliocentric coordinates
    Return: A possible i
    """

count = 0
i = np.array([180,0]) #i[0] is the minimum of i and i[1] is the maximum of i 
dRange = startObsDis(sNite=57048.068,sRa=85.621,sDec=11.499,eNite=57058.173,\
                     eRa=85.606,eDec=11.5013)
print(dRange)
#for d in dRange:
#    sP = startObsPoint(d,...)
#    eP = endObsPoint(sP)
#    iTemp = iFinding(sP, eP)
#    #print("d: ", d, "i: ", iTemp)
#    if iTemp < i[0]:
#        i[0] = iTemp
#    if iTemp > i[1]:
#        i[1] = iTemp
    
      

