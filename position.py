# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import math

"""
variables
"""
e = 0.4 # eccentricity

"""
Constants
"""
G = 6.67408e-11 # m^3/kg s gravitational constant
pi = np.pi
M = 1.989e30 #kg mass of the sun
k = 0.017202 # AU day M 0.01720209895 Gaussian gravitational constant

"""
global variables
"""
rad_array = np.zeros(361)

def period_calc(a):
    """
    Para: a (semi_axis length in AU)
    Return: tau (period in days)
    """
    massTotal = 1
    
    kep_k = 4*(pi**2)/((k**2)*(massTotal))
    
    tau = np.sqrt(kep_k*(a**3))
     
    return tau

def f_RHS(f):
    """
    evaluate the right hand side of the time-asteroid_angle equation
    """
    sin = math.sin(f)
    cos = math.cos(f)
    tan_half = math.tan(f/2)
    
    tan_outer = math.atan(tan_half * math.sqrt((1-e)/(1+e)))
    
    numerator = e*math.sqrt(1-e**2)*sin
    
    denominator  = 1 + e*cos
    
    result = 2*tan_outer - numerator/denominator
    
    return result

def f_fill():
    """
    filling the possible solution array with evaluations
    note possible error: wierod orientation of the ellipse
    """
    index = 0
   
    global rad_array
    
    for deg in np.linspace(0,2*pi,361):
        RHS = f_RHS(deg)
        if RHS < 0:
            RHS = 2*pi + RHS
        rad_array[index] = RHS
        print(index, deg, rad_array[index])
        index = index + 1 
        
    return
    
def f_calc(t, tau):
    """
    must call f_fill first!
    Para: t (current time in days)
    Return: f (asteroid angular position in rad)
    """
    
    LHS = 2*pi*t/tau
    
    root = 0
    diff = abs(LHS - rad_array[0])
    
    for index in np.arange(361):
        if(abs(LHS - rad_array[index]) < diff):
            root = index
            diff = abs(LHS - rad_array[index])
    
    return (root/360*2*pi)
    
f_fill()   
print(f_calc(0.79/2/pi,1))


    

