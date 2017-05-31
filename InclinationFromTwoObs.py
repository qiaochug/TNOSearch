#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 26 13:39:09 2017

@author: qiaochug
"""

import numpy as np
import matplotlib.pyplot as plt
import conversionsNoAst as cnv
from math import radians, degrees, atan2, asin, sin, cos, pi,sqrt,pi,acos
from vectorCalc import predict,maxVel,maxDegPerDay
diffs = 1
def startObsDis(sNite,sRa,sDec,eNite,eRa,eDec):
    """
    Para: start and end observation
    Return: range of probably distance
    """
    global diffs
    #divided the 20-600 AU into steps of distance
    steps = 40
    ds = np.linspace(1,600,steps)
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
    ranged = ([ds[max(0,radP-1)],ds[min(radP+1, steps -1)]])      
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

snite = 57277.063257
enite = 57278.026645
sra =  -50.5936171
era = -50.6109729
sdec = -46.843937
edec = -46.8387915
count = 0
numd = 20
i = np.array([180,0]) #i[0] is the minimum of i and i[1] is the maximum of i 
([mind,maxd]) = startObsDis(sNite=snite,sRa=sra,sDec=sdec,eNite=enite,\
                     eRa=era,eDec=edec)
dRange = np.linspace(mind,maxd,numd)
#print(dRange)
for d in dRange:
    sP = startObsPoint(sObsDis = d, sNite = snite, sRa = sra, sDec = sdec)
    #print("start Position: ",sP)
    eP = endObsPoint(startP = sP,eNite = enite,eRa = era,eDec = edec)
    #print("end Position: ", eP)
    iTemp = iFinding(startP = sP, endP = eP)
#    #print("d: ", d, "i: ", iTemp)
    if iTemp < i[0]:
        i[0] = iTemp
    if iTemp > i[1]:
        i[1] = iTemp
    print(iTemp)
print(i)
    
      

