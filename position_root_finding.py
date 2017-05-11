# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import math
import scipy

"""
variables
"""
a = 40 # semi-axis AU
e = 0 # eccentricity
i = 0 # inclination
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
phi = 0 # angle between the asteroid's position and the ascending node in the
#orbital plane
r = 0 #heliocentric distance of the asteroid
u = 0 # the projected heliocentric distance of an asteroid in the XY plane
alpha = 0 # angle between the ascending node and a line form the Sun to the 
#ecliptic plane above or below the asteroid's orbital position
Lambda_e = 0 #heliocentric longitude (kept stationary temporarily)
lam = 0 # lambda geometric ecliptic longitude
beta = 0 # geocentric ecliptic latitude
t = 0
tau = 0

def period_calc(a):
    """
    Para: a (semi_axis length in AU)
    Return: tau (period in days)
    """
    global tau
    
    massTotal = 1
    
    kep_k = 4*(pi**2)/((k**2)*(massTotal))
    
    tau = np.sqrt(kep_k*(a**3))
     
    return tau

def f_evaluate(f):
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
    
    RHS = 2*tan_outer - numerator/denominator
    
    if t%tau < tau/ 2:
        LHS = 2*pi*(t%tau)/tau
    
    else:
        LHS = 2*pi*(t%tau)/tau - 2*pi
    
    result = RHS - LHS
    
    return result

    
def f_calc():
    """
    must call f_fill first!
    Para: t (current time in days)
    Return: f (asteroid angular position in rad)
    """
#    if t%tau > tau/2:
#        #print("here")
#        root = scipy.optimize.newton(f_evaluate, 3/2*pi)
#    
#    if t%tau <= tau/2:
        #print("there")
    p = np.zeros(10)
    count = 0
    for num in np.linspace(0, 2*pi,10):
        try:
            root = scipy.optimize.newton(f_evaluate, num)
            p[count] = root
            count = count + 1
        except RuntimeError:
            pass
        
    #print(p)
    
    for num in np.arange(0,10):
        if p[num] >= 0 and p[num] <= 2*pi:
            root = p[num]
            
    
    #print(root)
    
    return (root)

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
    #print("r {:f}".format(r))
      
    u = r*math.sqrt((math.sin(phi))**2*(math.cos(i))**2+(math.cos(phi))**2)
    #print("u {:f}".format(u))
    
    alpha = math.atan(math.tan(phi)*math.cos(i))
    #print("alpha {:f}".format(alpha))
    
    Lambda_e = Lambda_e_calc(t)
    #print("lambda_e {:f}".format(Lambda_e))
    
    xyz[0] = u*math.cos(Omega + alpha) - R* math.cos(Lambda_e)
    #print("u {:f} cos(Lambda_e) {:f} cos(Omega + alpha) {:f}".format(u,math.cos(Lambda_e), math.cos(Omega + alpha)))
    
    xyz[1] = u*math.sin(Omega + alpha) - R* math.sin(Lambda_e)
    
    xyz[2] = r*math.sin(phi)*math.sin(i)
    #print("sin(phi) {:f} sin(i) {:f}".format(math.sin(phi), math.sin(i)))
    
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
    
    global t
    
    period_calc(a)
    
    days = np.arange(startDay, endDay+1)
    
#    largest = 0
#    largest_day = 0
#    negative = False
#    f_large = False
#    f_negative = False
    
    for day in days:
        t = day
        f = f_calc()
        #print("day {:f} tau {:f}".format(day, tau))
        xyz = xyz_calc(f, day)
        #print("x {:f} y {:f} z {:f}".format(xyz[0], xyz[1], xyz[2]))
        lb = lambda_beta_calc(xyz)
#        if lb[0] < 0:
#            lb[0] = lb[0] + 2*pi
        print("day {:d} f {:f} lambda {:f} Beta {:f}".format(day, f,lb[0]/2/pi*360, lb[1]/2/pi*360))
#        if lb[1]< largest:
#            largest = lb[1]
#            largest_day = day
#        if lb[1] < 0:
#            negative = True
#        if f < 0:
#            f_negative = True
#        if f > 4:
#            if f < 7:
#                f_large = True
            
    #print(largest_day, largest, negative, f_large, f_negative)  
    return
 
generate(0,365)       
#f_fill()   
#print(f_calc(0.6,1))
#xyz = xyz_calc(2,30)
#print(xyz[0], xyz[1], xyz[2])
#lb = lambda_beta_calc(xyz)
#print(lb[0],lb[1])



    

