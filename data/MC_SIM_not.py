# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 16:58:57 2022

@author: jbseo
"""

import numpy as np
import matplotlib.pyplot as plt



def getAcc( vel, pos ):
    theta = np.arctan2(vel[:,1],vel[:,0])

    vv = np.sqrt(vel[:,1]**2+vel[:,0]**2)
    
    a = np.zeros((N,2))
    a[:,0] = vv*B*np.cos(theta+0.5*np.pi)
    a[:,1] = vv*B*np.sin(theta+0.5*np.pi)
    return a


def bound( pos ):
    pos[:,0][pos[:,0]>1] -= 1
    pos[:,0][pos[:,0]<0] += 1
    pos[:,1][pos[:,1]>1] -= 1
    pos[:,1][pos[:,1]<0] += 1
    return pos

# Simulation parameters
N         = 10000    # Number of particles
t         = 0      # current time of the simulation
tEnd      = 1   # time at which simulation ends
dt        = 0.01   # timestep
B         = 50

# Generate Initial Conditions
np.random.seed()            # set the random number generator seed

pos  = np.random.rand(N,2)   # randomly selected positions and velocities
vel  = np.ones((N,2))
vel[:,0]  = np.sin(pos[:,1]*np.pi*4)
vel[:,1]  = 0

# calculate initial gravitational accelerations
acc = getAcc( vel, pos )


# number of timesteps
Nt = int(np.ceil(tEnd/dt))

# prep figure


for i in range(Nt):
    vel += acc * dt/2.0 # (1/2) kick
    pos += vel * dt # drift
    acc = getAcc( vel, pos ) # update accelerations
    vel += acc * dt/2.0 # (1/2) kick
    t += dt # update time
    pos = bound(pos)
    
    fig = plt.figure(figsize=(5,5), dpi=80)
    plt.scatter(pos[:,0],pos[:,1],s=10,color='blue')
    plt.xlim([0, 1]); plt.ylim([0, 1])
    plt.savefig('data%03d.png'%i)
    plt.show()
        
