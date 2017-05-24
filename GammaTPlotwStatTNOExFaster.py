#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 13:41:05 2017

@author: qiaochug
"""

import numpy as np
import matplotlib.pyplot as plt
from math import atan2,cos,sin,sqrt

import matplotlib.pyplot as plt
from matplotlib import animation


pi = np.pi
d = 20 # d is the hypothetical distance from earth to TNO at the starting
        # observation, it is set to such that the TNO has the greatest apparent
        # motion in the sky(to exaggerate the motion, we could further reduce the distance)

def predict(gammaS, betaS, deltaT):
    """
    Para: Gamma and beta at starting date, and time difference(in days)
    Return: New Gamma and New Beta(could be new lambda and new beta,
    see ******************************)
    Note: the Lambda Earth is approximated using the period of Earth, 
    the distance is in units of AU, the coordinates system is built with
    sun as the center, earth starts at (0,-1) and ecliptic plane as the 
    X-Y plane
    Do not omit the global variables above!
    """
    
    #Lambda Earth(in fact, the difference in lambda earth) 
    #is "approximated" using the Earth Period
    lambdaE = (deltaT/365.25)*2*pi
    
    #result[0] = new gamma, result[1] = new beta
    result = np.array([0.0,0.0])
    
    dp = d*cos(betaS)
    Xt = dp*sin(gammaS)
    Yt = -1-dp*cos(gammaS)
    Zt = d*sin(betaS)
    Xe = sin(lambdaE)
    Ye = -cos(lambdaE)
    #Ze = 0
    
    # calculation for new Gamma
    # vector(Xt-Xe,Yt-Ye)
    gammaN = atan2(Yt-Ye, Xt-Xe) - lambdaE + pi/2
    #gammaN = atan2(Yt-Ye, Xt-Xe) + pi/2 #*************************For beta against lambda just change 
    #the line above to this code**************************
    if gammaN < 0:
        gammaN = gammaN + 4*pi
    gammaN = gammaN %(2*pi)
    
    result[0] = gammaN
    #calculation for new Beta
    proj = sqrt((Xt-Xe)**2 + (Yt-Ye)**2)
    
    result[1] = atan2(Zt, proj)
    
    return result

######################### Everything below are for testing purpose#######################
def to_rad(deg):
    return deg/360*2*pi

def to_deg(rad):
    return rad/(2*pi)*360

count = 0
x_steps = np.zeros(360)#gamma
y_steps = np.zeros(360)#beta
times = np.linspace(0,365.25,360)
for t in times:
    re = predict(gammaS = 5.1065, betaS = -0.4788, deltaT = t)
    print(re)
    x_steps[count] = to_deg(re[0]) #+ (t/365.25)*2*pi)%(2*pi)
    y_steps[count] = to_deg(re[1])
    count = count + 1

#plt.close('all')
#fig = plt.figure()
#ax = plt.axes(xlim=(30,50), ylim=(10,30))
#
#(my_line,) = ax.plot([],[],lw = 2)
#(my_point,) = ax.plot([],[],'ro',ms = 3)
#
#def get_step(n,x,y,this_line,this_point):
#    this_line.set_data(x[:n+1], y[:n+1])
#    this_point.set_data(x[n], y[n])
#
#
#my_movie = animation.FuncAnimation(fig,get_step, frames= 2,\
#                                   fargs= (x_steps,y_steps,my_line,my_point))
#my_movie.save('GammaAgnTEx5.mp4',fps = 5)