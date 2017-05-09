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
a = 40 # semi-axis AU
e = 0 # eccentricity
i = 0.5 # inclination
Omega = 0 # longitude of the ascending node
omega = 0 #argument of perihelion

"""
Constants
"""
G = 6.67408e-11 # m^3/kg s gravitational constant
pi = np.pi
M = 1.989e30 #kg mass of the sun
k = 0.017202 # AU day M 0.01720209895 Gaussian gravitational constant
R = 1 # AU earth-sun distance

"""
global variables
"""
rad_array = np.zeros(361)
phi = 0 # angle between the asteroid's position and the ascending node in the
#orbital plane
r = 0 #heliocentric distance of the asteroid
u = 0 # the projected heliocentric distance of an asteroid in the XY plane
alpha = 0 # angle between the ascending node and a line form the Sun to the 
#ecliptic plane above or below the asteroid's orbital position
Lambda_e = 0 #heliocentric longitude (kept stationary temporarily)
lam = 0 # lambda geometric ecliptic longitude
beta = 0 # geocentric ecliptic latitude

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
    helper function, called in f_fill()
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
    helper function, user should call it before the first f_calc
    filling the possible solution array with evaluations
    note possible error: negative evaluations at angles larger than pi
    """
    index = 0
   
    global rad_array
    
    for deg in np.linspace(0,2*pi,361):
        RHS = f_RHS(deg)
        if RHS < 0:
            RHS = 2*pi + RHS
        rad_array[index] = RHS
        #print(index, deg, rad_array[index])
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

def Lambda_e_calc(t):
    """
    helper function, called in xyz_calc
    update earth heliocentric longitude
    """
    Le = t/365.25 * 2 *pi
    return Le

def xyz_calc(f,t):
    """
    compute xyz coordinates of the asteroid
    Return: an array of size 3, x y z
    """
    
    xyz = np.array([0.0,0.0,0.0])
    
    global phi
    global r
    global u
    global alpha
    global Lambda_e
    
    phi = omega + f
    
    r = a*(1-e**2)/(1+e*math.cos(f))
    
    u = r*math.sqrt((math.sin(phi))**2*(math.cos(i))**2+(math.cos(phi))**2)
    #print("u {:f}".format(u))
    
    alpha = math.atan(math.tan(phi)*math.cos(i))
    #print("alpha {:f}".format(alpha))
    
    Lambda_e = Lambda_e_calc(t)
    #print("lambda_e {:f}".format(Lambda_e))
    
    xyz[0] = u*math.cos(Omega + alpha) - R* math.cos(Lambda_e)
    
    xyz[1] = u*math.sin(Omega + alpha) - R* math.sin(Lambda_e)
    
    xyz[2] = r*math.sin(phi)*math.sin(i)
    
    #print("x {:f} y {:f} z {:f}".format(xyz[0], xyz[1], xyz[2]))
    
    return xyz

def lambda_beta_calc(xyz):
    
    lb = np.array([0.0,0.0])
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]
    #print("x {:f} y {:f} z {:f}".format(x, y, z))
    
    lb[0] = math.atan(y/x)
    #print(math.atan(y/x))
    #print("lb[0] {:f}".format(lb[0]))
    lb[1] = math.atan(z/math.sqrt(x**2+y**2))
    #print("l {:f} b {:f}".format(lb[0],lb[1]))
    
    return lb


def generate(startDay, endDay):
    """
    print out the lambda beta trace given 
    a period of time
    """
    
    tau = period_calc(a)
    
    f_fill()
    
    days = np.arange(startDay, endDay+1)
    
    for day in days:
        f = f_calc(day, tau)
        #print("day {:f} tau {:f}".format(day, tau))
        xyz = xyz_calc(f, day)
        #print("x {:f} y {:f} z {:f}".format(xyz[0], xyz[1], xyz[2]))
        lb = lambda_beta_calc(xyz)
        #if lb[0] < 0:
            #lb[0] = lb[0] + 2*pi
        print("day {:d} f {:f} lambda {:f} Beta {:f}".format(day, f,lb[0]/2/pi*360, lb[1]/2/pi*360))
        
    return
 
generate(0,365)       
#f_fill()   
#print(f_calc(0.6,1))
#xyz = xyz_calc(2,30)
#print(xyz[0], xyz[1], xyz[2])
#lb = lambda_beta_calc(xyz)
#print(lb[0],lb[1])



    

