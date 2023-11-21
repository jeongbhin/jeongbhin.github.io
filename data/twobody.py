# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 17:01:18 2021

@author: jbseo
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt


def a_cal( pos, mass ):
    a = np.zeros((2,2))
    for i in range(2):
        for j in range(2):
            if i != j :
                dx = pos[j,0] - pos[i,0]
                dy = pos[j,1] - pos[i,1]
                r = (dx**2 + dy**2 +rc**2)**(0.5)
                a[i,0] +=  G * (dx * r**(-3)) * mass[j]
                a[i,1] +=  G * (dy * r**(-3)) * mass[j]
    return a

def dec( pos, vel,dt ):
    dx = pos[1,0] - pos[0,0]
    dy = pos[1,1] - pos[0,1]
    r = (dx**2 + dy**2)**(0.5)
    if r <= 3/40*0.5 :
        vel = vel*(1-dt/t_dcay)
    return vel

t         = 0     
t_end      = 0.5e0 
dt        = 3e-3
G         = 1.05 
rc = 3/40
t_dcay = 0.03
mass = np.array([0.5,1.0])
mass = np.reshape(mass,(2,1))
pos  = np.array([[-0.2,-rc],[0.1,rc]])
vel  = np.array([[1,-0.0],[-0.0,0.0]])

#vel -= np.mean(mass * vel) / np.mean(mass)
acc = a_cal( pos, mass )
		
Nt = int(t_end/dt)	


for i in range(Nt):    
    
    vel += acc * dt/2.0
    
    pos += vel * dt
    acc = a_cal( pos, mass ) 
    vel += acc * dt/2.0
    vel = dec(pos,vel,dt)
    t += dt
    
    if i%3 == 0:
        fig = plt.figure(figsize=(5,5))
        
        plt.scatter(pos[:,0],pos[:,1],s=mass*100,color='blue',alpha=0.5)
        plt.xlim([-0.5,0.5])
        plt.ylim([-0.5,0.5])
        plt.savefig('%03d.png'%i)
        plt.show()