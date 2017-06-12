#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 26 13:39:09 2017

@author: qiaochug
"""
import time
import sys
sys.path.append('/global/project/projectdirs/dessn/diffim/edison/destnos/pyOrbfit/') 
import ephem 
import pickle
from astropy.time import Time  
from LinkerLib import Detection
from Orbit import Orbit
from PartitionSeason import printPercentage
from LinkerLib import Triplet
import argparse
import numpy as np
import conversionsNoAst as cnv
from math import radians, degrees, atan2, asin, sin, cos, sqrt,pi,acos
from vectorCalc import predict,maxVel,maxDegPerDay
from datetime import datetime

diffs = 1
def startObsDis(sNite,sRa,sDec,eNite,eRa,eDec):
    """
    Para: start and end observation
    Return: range of probably distance
    """
    global diffs
    #divided the 20-600 AU into steps of distance
    steps = 40
    ds = np.linspace(15,600,steps)
    diffs = np.zeros(steps)
    Index = 0
    for d in ds:
        #compute the vector
        (betaS, gammaS) = cnv.equatorialToEcliptic(sRa, sDec, sNite)
        gammaN, betaN = predict(np.radians(gammaS), np.radians(betaS), (eNite - sNite),d)
        (raN, decN) = cnv.eclipticToEquatorial(degrees(betaN), degrees(gammaN), eNite)
        if raN > 180:
            raN = raN - 360
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
        radP = int(round((firstZInd + lastZInd)/2))
    # take 3 values around the radiation points as range of d
    # consider the edge cases where the radP is smaller than 2 or larger than steps - 2
    #print(radP)
    ranged = ([ds[max(0,radP-1)],ds[min(radP+1, steps-1)]])      
    return ranged
    
def startObsPoint(sObsDis, sNite, sRa,sDec):
    """
    Para: start observation, and suspected distance(one)
    Return: A point in heliocentric coordinates
    """
    (betaS, gammaS) = cnv.equatorialToEcliptic(sRa, sDec, sNite)
    lambdaE = ((sNite - 56192.447049) / 365.25 * 360) % (360)
    lam = (lambdaE + gammaS) % (360)
    lam = radians(lam)
    betaS = radians(betaS)
    lambdaE = radians(lambdaE)
    #projected distance
    pd = sObsDis * cos(betaS)
    #X Y Z in geocentric coordinates
    Xe = pd * cos(lam)
    Ye = pd * sin(lam)
    Ze = sObsDis * sin(betaS)
    Xs = Xe + cos(lambdaE)
    Ys = Ye + sin(lambdaE)
    return ([Xs,Ys,Ze])
    
def endObsPoint(startP, eNite,eRa,eDec):
    """
    Para: start point in heliocentric coordinates, end observation
    Return: A point in heliocentric coordinates
    """
    sX = startP[0]
    sY = startP[1]
    sZ = startP[2]
    #position of the earth in heliocentric coordinates
    lambdaE = ((eNite - 56192.447049) / 365.25 * 360) % (360)
    lambdaE = radians(lambdaE)
    Ex = cos(lambdaE)
    Ey = sin(lambdaE)
    d_hypo = sqrt((sX-Ex)**2+(sY-Ey)**2+(sZ)**2)
    
    #end night lambda and beta
    (betaE, gammaE) = cnv.equatorialToEcliptic(eRa, eDec, eNite)
    betaE = radians(betaE)
    gammaE = radians(gammaE)
    lam = (lambdaE + gammaE) % (2*pi)
    #projected distance
    pd = d_hypo * cos(betaE)
    # X Y Z in geocentric coordinates
    Xe = pd * cos(lam)
    Ye = pd * sin(lam)
    Ze = d_hypo * sin(betaE)
    # XYZ in heliocentric coordinates
    Xs = Xe + cos(lambdaE)
    Ys = Ye + sin(lambdaE)
    return ([Xs,Ys,Ze])
    
def iFinding(startP, endP):
    """
    Para: start point and end point in heliocentric coordinates
    Return: A possible i
    """
    xs = startP[0]
    ys = startP[1]
    zs = startP[2]
    
    xe = endP[0]
    ye = endP[1]
    ze = endP[2]
    
    # the normal vector(Vs X Ve), order matters
    xn = ys*ze-zs*ye
    yn = zs*xe-xs*ze
    zn = xs*ye-ys*xe
    
    # unit vector normal to the ecliptic plane
    xu = 0
    yu = 0
    zu = 1
    
    # dot product of two normal vectors
    dotProd = xn*xu + yn*yu + zn*zu
    #||x||
    lenOfn1 = sqrt((xn)**2+(yn)**2+(zn)**2)
    lenOfn2 = sqrt((xu)**2+(yu)**2+(zu)**2)
    inc = acos(dotProd/(lenOfn1 * lenOfn2))
    return degrees(inc)

def calcInc(pairs, outputFile):
    with open(outputFile, 'w+') as output:
        size = len(pairs)
        counter = size 
        starT = time.time()
        # for each detection, loop through each possible pairs
        for det in pairs:
            printPercentage(size - counter, size, time.time()-starT)
            for link in det.linkedList:
                #make the triplet object
                triplet = Triplet([det, link])
                orbit = triplet.calcOrbit() 
                elements, errs = orbit.get_elements()
                #output.write('raLst ' + str(raLst) + '\n')
                #output.write('decLst ' + str(decLst) + '\n') 
                #output.write('dateLst ' + str(dateLst) + '\n') 
                if elements['a'] >= 20:
                    output.write('det1: ' + det.toStr() + '\n')
                    output.write('det2: ' + link.toStr() + '\n')  
                    #output.write('elements: ' + str(elements) + '\n')
                    #output.write('errs: ' + str(errs) + '\n')
                    
                    #my code of inclination
                    snite = det.mjd
                    enite = link.mjd
                    sra =  det.ra
                    era = link.ra
                    sdec = det.dec
                    edec = link.dec
                    numd = 20
                    i = np.array([180.0,0.0]) #i[0] is the minimum of i and i[1] is the maximum of i 
                    ([mind,maxd]) = startObsDis(sNite=snite,sRa=sra,sDec=sdec,eNite=enite,\
                                         eRa=era,eDec=edec)
                    dRange = np.linspace(mind,maxd,numd)
                    #print(dRange)
                    tot = 0
                    for d in dRange:
                        sP = startObsPoint(sObsDis = d, sNite = snite, sRa = sra, sDec = sdec)
                        #print("start Position: ",sP)
                        eP = endObsPoint(startP = sP,eNite = enite,eRa = era,eDec = edec)
                        #print("end Position: ", eP)
                        iTemp = iFinding(startP = sP, endP = eP)
                        tot = tot + iTemp
                    #    #print("d: ", d, "i: ", iTemp)
                        if iTemp < i[0]:
                            i[0] = iTemp
                        if iTemp > i[1]:
                            i[1] = iTemp
                        #print(iTemp)
                    #output.write(str(snite)+'' +str(enite)+' '+ str(sra)+' '+str(era)+' '+str(sdec)\
                     #            +' '+str(edec)
                    mean = tot/numd
                    if abs(i[0] - mean) > abs(i[1]-mean):
                        smarterPred = (i[1]+mean)/2
                    else:
                        smarterPred = (i[0]+mean)/2
                    output.write('OrbitFit A: '+ str(elements['a'])+'\n')
                    output.write('Planar Method Estimation of d: '+ str(dRange[0])+ ' to '\
                                                                        + str(dRange[numd-1])+'\n')
                    output.write('OrbitFit I: '+ str(elements['i'])+' OrbitFit Err: ' + str(errs['i'])+'\n')
                    output.write('Planar Method Estimation of i: '+ str(i)+' Smartmean: '+str(smarterPred)+'\n')
                    output.write('*****************************\n\n\n')
            counter -= 1
             
def main():
    args = argparse.ArgumentParser()
    args.add_argument('linkedPairs', nargs=1, 
                        help='path to pickle file')
    args.add_argument('outFile', nargs=1,
                        help='path to output')
    args = args.parse_args()
    
    #load the pickle file with list of detections and links
    detPairs = pickle.load(open(args.linkedPairs[0], 'rb'))
    #calculate the orbits of potential triplets
    calcInc(detPairs, args.outFile[0])
    
if __name__ == '__main__':
    main()
 
      

