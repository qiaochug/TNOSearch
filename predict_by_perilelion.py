#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 17:06:26 2017

@author: qiaochug
"""

import numpy as np
import math
import scipy
import matplotlib.pyplot as plt
from matplotlib import animation
import astropy.time
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, EarthLocation
from astropy.coordinates import get_body_barycentric, get_body
import astropy.coordinates
loc = EarthLocation.of_site('greenwich')

"""
variables
"""
a = 40 # semi-axis AU
e = 0 # eccentricity
i = 0.4 # inclination
Omega = 0 # longitude of the ascending node
omega = 0 #argument of perihelion

start_f = 0 # the starting offset for asteroid
start_date = Time("2000-12-30")
perihelion_date = Time("2000-12-30")
date = Time("2000-12-30")

"""
Constant
"""
G = 6.67408e-11 # m^3/kg s gravitational constant
pi = np.pi
M = 1.989e30 #kg mass of the sun
k = 0.017202 # AU day M 0.01720209895 Gaussian gravitational constant
R = 1 # AU earth-sun distance
tilt = 0.4088

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
    helper function, called by f_calc()
    evaluate the difference between right hand side 
    and left hand side of the time-asteroid_angle equation
    for minimization
    """
    sin = math.sin(f)
    cos = math.cos(f)
    tan_half = math.tan(f/2)
    
    tan_outer = math.atan(tan_half * math.sqrt((1-e)/(1+e)))
    
    numerator = e*math.sqrt(1-e**2)*sin
    
    denominator  = 1 + e*cos
    
    RHS = 2*tan_outer - numerator/denominator
    
    LHS = 2*pi*(t%tau)/tau
    
    if RHS < 0:
        RHS = RHS + 2*pi
        
    result = RHS - LHS
    
    return abs(result)
    
    
def f_calc():
    """
    t (current time in days) accessed through global variales
    Return: f (asteroid angular position in rad)
    """
    p = np.zeros(10)
    count = 0
    for num in np.linspace(0, 2*pi,10):
        try:
            root = scipy.optimize.newton(f_evaluate, num)
            p[count] = root
            count = count + 1
        except RuntimeError:
            pass
    
    for num in np.arange(0,10):
        if p[num] > 0 and p[num] <= 2*pi:
            root = p[num]
            
    root = (root + start_f)%(2*pi)
    
    return (root)

def Lambda_e_calc(date):
    """
    helper function, called in xyz_calc
    update earth heliocentric longitude
    using astropy library data
    """
    date = Time(date, format = 'jd')
    
    with solar_system_ephemeris.set('builtin'):
        earth = get_body('earth',date,loc)
    
    earthAngle = earth.transform_to('heliocentrictrueecliptic').lon.value

    return earthAngle

def xyz_calc(f,date):
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
    
    phi = (omega + f)%(2*pi)
    
    r = a*(1-e**2)/(1+e*math.cos(f))
      
    u = r*math.sqrt(((math.sin(phi))**2)*((math.cos(i))**2)+(math.cos(phi))**2)
    
    to_be_tan = math.tan(phi)*math.cos(i)
    tan_1 = math.atan2(to_be_tan, 1)
    if tan_1 < 0:
            tan_1 = tan_1 + 2*pi
    tan_2 = math.atan2(-to_be_tan, -1)
    if tan_2 < 0:
            tan_2 = tan_2 + 2*pi
    
    if phi >=0 and phi <= pi:
        if tan_1 >= 0 and tan_1 <= pi:
            alpha = tan_1
        else:
            alpha = tan_2
    
    if phi > pi and phi <= 2*pi:
        if tan_1 > pi and tan_1 <= 2*pi:
            alpha = tan_1
        else:
            alpha = tan_2
            
    Lambda_e = Lambda_e_calc(date)
    
    xyz[0] = u*math.cos(Omega + alpha) - R* math.cos(Lambda_e)
    
    xyz[1] = u*math.sin(Omega + alpha) - R* math.sin(Lambda_e)
    
    xyz[2] = r*math.sin(phi)*math.sin(i)
    
    return xyz

def lambda_beta_calc(xyz):
    """
    legacy code for lamda and beta calculation, not used in predict_peri
    lamda and beta are the position of asteroid in respect to earth
    parameter: [x,y,z]
    return: [lamda,beta]
    """
    lb = np.array([0.0,0.0])
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]
    
    lb[0] = math.atan2(y,x)
    if lb[0] < 0:
        lb[0] = lb[0] + 2*pi                     

    lb[1] = math.atan(z/math.sqrt(x**2+y**2))
    
    return lb

def to_rad(deg):
    return deg/360*2*pi

def to_deg(rad):
    return rad/(2*pi)*360

def parse_date(string):
    """
    parse in the date which has the format "YYYY-MM-DD.dddd"
    to "YYYY-MM-DD HH:MM:SS:ssss" for Time object initialization
    """
    
    DT = string.split(".")
    day = DT[0]
    time = (float)("0." +str(DT[1]))
    total_seconds = time * 24*60*60
    hour = (int)(total_seconds/(60*60))
    minute = (int)((total_seconds- hour*3600)/60)
    second = total_seconds-hour*3600-minute*60
    if hour < 10:
        hour = "0"+ str(hour)
    else:
        hour = str(hour)
    
    if minute < 10:
        minute = "0" + str(minute)
    else:
        minute = str(minute)
    
    if second < 10:
        second = "0" + str(second)
    else:
        second = str(second)
    
    new_string = day + " " + hour+":"+minute+":"+second
    return new_string  

#enter i Omega and omega in degrees
#start date and perihelion date in format "YYYY-MM-DD.dddd"
#reqeust date in format ["YYYY-MM-DD","YYYY-MM-DD"...]
def predict_peri(a_p,e_p,i_p,Omega_p,omega_p,start_date_p,perihelion_date_p, request_dates):
    """
    predict the position of the asteroid given a,e,i,Omega,omega 
    and perihelion date (pd used to calculate true anamoly)
    display position in ra(HH MM SS) and dec(DD MM SS)
    """
    global a 
    global e
    global i
    global Omega
    global omega
    global start_date
    global perihelion_date
    global date
    global t
    
    #initialize variables
    a = a_p
    e = e_p
    i = to_rad(i_p)
    Omega = to_rad(Omega_p)
    omega = to_rad(omega_p)
    period_calc(a)
    perihelion_date = Time(parse_date(perihelion_date_p))
    
    #compute start date true anomaly
    start_date = Time(parse_date(start_date_p))
    date = start_date
    t = perihelion_date.jd - date.jd
    t = tau - t
    f = f_calc()
    
    #compute start date position coordinates
    xyz = xyz_calc(f, date)
    coords = astropy.coordinates.CartesianRepresentation(x=xyz[0],y = xyz[1],z = xyz[2], unit = 'AU')
    new_coord = astropy.coordinates.SkyCoord(coords, frame = 'heliocentrictrueecliptic')
    new_coord = new_coord.transform_to('gcrs')
    
    #convert to display format
    ra = new_coord.ra.value
    hour = (int)(ra/360*24)
    minute = (int)((ra - 15*hour)/15*60)
    second = (int) ((ra - 15*hour-0.25*minute)/0.25*60)
    print(start_date_p,hour,"h",minute,"m",second,"s", new_coord.dec)
    
    #iterate throught the date list and compute, display the coordinates of asteroid on each date
    for date_itr in request_dates:
        date = Time(parse_date(date_itr))
        t = perihelion_date.jd - date.jd
        t = tau - t
        
        f = f_calc()
        xyz = xyz_calc(f, date)
        coords = astropy.coordinates.CartesianRepresentation(x=xyz[0],y = xyz[1],z = xyz[2], unit = 'AU')
        new_coord = astropy.coordinates.SkyCoord(coords, frame = 'heliocentrictrueecliptic')
        new_coord = new_coord.transform_to('gcrs')
        
        ra = new_coord.ra.value
        hour = (int)(ra/360*24)
        minute = (int)((ra - 15*hour)/15*60)
        second = (int) ((ra - 15*hour-0.25*minute)/0.25*60)
        print(date_itr,hour,"h",minute,"m",second,"s", new_coord.dec)
        
        
    return

#user only need to use the predict_peri function, by modifying the parameters below
days = np.array(["2014-09-27.311201","2014-10-21.357384","2016-07-18.413576"])
predict_peri(a_p = 108.7787824,e_p = 0.6514653,i_p = 26.78539,Omega_p = 131.08649,omega_p = 29.68223,\
             start_date_p = "2014-08-19.424583",\
        perihelion_date_p = "2142-01-02.38251",request_dates = days)

