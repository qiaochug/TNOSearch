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

"""
variables
"""
a = 40 # semi-axis AU
e = 0 # eccentricity
i = 0.4 # inclination
Omega = 0 # longitude of the ascending node
omega = 0 #argument of perihelion

start_f = 0 # the starting offset for asteroid

"""
Constant
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
    
    LHS = 2*pi*(t%tau)/tau
    #print("LHS {:f}".format(LHS))
    
    if RHS < 0:
        RHS = RHS + 2*pi
    
    #print("RHS {:f}".format(RHS))
        
    result = RHS - LHS
    
    return abs(result)

#period_calc(a)
#for testangle in np.linspace(0,2*pi, 30):
#    print("angle {:f}".format(testangle))
#    f_evaluate(testangle)
    
    
def f_calc():
    """
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
    
    for num in np.arange(0,10):
        if p[num] > 0 and p[num] <= 2*pi:
            root = p[num]
            
    #print("t {:d} root {:f} p {}".format(t, root,p))
    
    root = (root + start_f)%(2*pi)
    
    return (root)

def Lambda_e_calc(t):
    """
    helper function, called in xyz_calc
    update earth heliocentric longitude
    """
    Le = (t%365.25)/365.25 * 2 *pi
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
    
    phi = (omega + f)%(2*pi)
    
    r = a*(1-e**2)/(1+e*math.cos(f))
    #print("r {:f}".format(r))
      
    u = r*math.sqrt(((math.sin(phi))**2)*((math.cos(i))**2)+(math.cos(phi))**2)
    #print("u {:f}".format(u))
    
    to_be_tan = math.tan(phi)*math.cos(i)
    #print(to_be_tan)
    tan_1 = math.atan2(to_be_tan, 1)
    #print(tan_1)
    if tan_1 < 0:
            tan_1 = tan_1 + 2*pi
    tan_2 = math.atan2(-to_be_tan, -1)
    #print(tan_2)
    if tan_2 < 0:
            tan_2 = tan_2 + 2*pi
            
    #print("tan_1 {:f} tan_2 {:f}".format(tan_1, tan_2))
   
   
    
    if phi >=0 and phi <= pi:
        if tan_1 >= 0 and tan_1 <= pi:
            #print("above1")
            alpha = tan_1
        else:
            #print("above2")
            alpha = tan_2
    
    if phi > pi and phi <= 2*pi:
        if tan_1 > pi and tan_1 <= 2*pi:
            #print("below1 {:f} {:f}".format(tan_1,tan_2))
            alpha = tan_1
        else:
            #print("below2")
            alpha = tan_2
    
    #alpha = math.atan(to_be_tan)
    #print("alpha {:f}".format(alpha))
    
    Lambda_e = Lambda_e_calc(t)
    #print("t {:d} lambda_e {:f}".format(t,Lambda_e))
    
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
    
    lb[0] = math.atan2(y,x)
    if lb[0] < 0:
        lb[0] = lb[0] + 2*pi
    #lb[0] = 2*pi - lb[0]                      

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
    
    largest = 0
    largest_day = 0
    negative = False
    f_large = False
    f_negative = False
    
    for day in days:
        t = day
        f = f_calc()
        #print("f {:f} day {:f}".format(f,day))
        #print("day {:f} tau {:f}".format(day, tau))
        xyz = xyz_calc(f, day)
        #print("day {:d} x {:f} y {:f} z {:f}".format(day, xyz[0], xyz[1], xyz[2]))
        lb = lambda_beta_calc(xyz)
#        if lb[0] < 0:
#            lb[0] = lb[0] + 2*pi
        
        #print("day {:d} f {:f} lambda {:f} Beta {:f}".format(day, f,lb[0]/2/pi*360, lb[1]/2/pi*360))
        if lb[1]< largest:
            largest = lb[1]
            largest_day = day
        if lb[1] < 0:
            negative = True
        if f < 0:
            f_negative = True
        if f > 4:
            if f < 7:
                f_large = True
            
    #print(largest_day, largest, negative, f_large, f_negative)  
    return

#generate(180,185)      
def animate(startDay, endDay):
    
    plt.close('all')
    bound = 365
    fig = plt.figure()
    ax = plt.axes(xlim=(-bound,bound), ylim=(-bound, bound))
    (my_line,) = ax.plot([],[],lw = 2)
    (my_point,) = ax.plot([],[],'ro',ms = 9)
    
    global t
    
    period_calc(a)
    
    days = np.arange(startDay, endDay+1)
    
    x_steps = np.zeros(endDay - startDay +1)
    y_steps = np.zeros(endDay - startDay +1)
        
    count = 0
    
    for day in days:
        t = day
        f = f_calc()
        #print("day {:f} tau {:f}".format(day, tau))
        xyz = xyz_calc(f, day)
        #print("x {:f} y {:f} z {:f}".format(xyz[0], xyz[1], xyz[2]))
        x = xyz[0]
        y = xyz[1]
        z = xyz[2]
        lb = lambda_beta_calc(xyz)
#        x_steps[count] = math.sin(lb[0]) * (math.sqrt(x**2+y**2))*10
#        y_steps[count] = math.sin(lb[1]) * (math.sqrt(x**2+y**2+z**2))*10
        x_steps[count] = count
        y_steps[count] = count
        count = count + 1
    
    print(x_steps)
    print(y_steps)
    
    def get_step(n,x,y,this_line,this_point):
        this_line.set_data(x[:n+1], y[:n+1])
        this_point.set_data(x[n], y[n])
    
    my_movie = animation.FuncAnimation(fig,get_step, frames=endDay - startDay +1,\
                                       fargs= (x_steps,y_steps,my_line,my_point))


plt.close('all')
fig = plt.figure()
ax = plt.axes(xlim=(0,5), ylim=(-3,3))

(my_line,) = ax.plot([],[],lw = 2)
(my_point,) = ax.plot([],[],'ro',ms = 3)

period_calc(a)

days = np.arange(1,1000)

x_steps = np.zeros(1000)
y_steps = np.zeros(1000)

count  = 0
for day in days:
        t = day
        f = f_calc()
        #print("day {:f} tau {:f}".format(day, tau))
        xyz = xyz_calc(f, day)
        #print("x {:f} y {:f} z {:f}".format(xyz[0], xyz[1], xyz[2]))
        x = xyz[0]
        y = xyz[1]
        z = xyz[2]
        lb = lambda_beta_calc(xyz)
#        x_steps[count] = math.sin(lb[0])* (math.sqrt(x**2+y**2))*20
#        y_steps[count] = math.sin(lb[1])* (math.sqrt(x**2+y**2+z**2))*20
        x_steps[count] = lb[0]/2/pi*360
        y_steps[count] = lb[1]/2/pi*360
#        x_steps[count] = y
#        y_steps[count] = x
        #x_steps[count] = count
        #y_steps[count] = count
        count = count + 1

print(x_steps)
print(y_steps)

def get_step(n,x,y,this_line,this_point):
    this_line.set_data(x[:n+1], y[:n+1])
    this_point.set_data(x[n], y[n])


my_movie = animation.FuncAnimation(fig,get_step, frames=999,\
                                   fargs= (x_steps,y_steps,my_line,my_point))
#my_movie.save('newphiAtane0.4omega0.mp4',fps = 30)

#generate(0,365)
#f_fill()   
#print(f_calc(0.6,1))
#xyz = xyz_calc(2,30)
#print(xyz[0], xyz[1], xyz[2])
#lb = lambda_beta_calc(xyz)
#print(lb[0],lb[1])



    

