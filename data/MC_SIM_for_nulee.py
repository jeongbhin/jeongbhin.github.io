# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 16:58:57 2022

@author: jbseo
"""

import numpy as np
import matplotlib.pyplot as plt


def getAcc( vel, pos ):
    theta = np.arctan2(vel[:,1],vel[:,0])

    vv = np.sqrt(vel_ini[:,1]**2+vel_ini[:,0]**2)
    
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


def histogram(pos):
    hist,bins = np.histogram(pos[:,1],bins=50)
    bins = 0.5*(bins[1:]+bins[:-1])

    return hist,bins

# Simulation parameters
N         = 10000    # Number of particles
t         = 0      # current time of the simulation
m         = 1
q         = 1
B         = 100

dt        = m/(q*B)*0.01   # timestep


tEnd      = 1   # time at which simulation ends

# Generate Initial Conditions
np.random.seed()            # set the random number generator seed

pos  = np.random.rand(N,2)   # randomly selected positions and velocities
vel  = np.ones((N,2))
vel[:,0]  = np.sin(pos[:,1]*np.pi*4) 
vel[:,1]  = 0
vel_ini = vel

# calculate initial gravitational accelerations
acc = getAcc( vel, pos )


# number of timesteps
Nt = int(np.ceil(tEnd/dt))

fig = plt.figure(figsize=(20,5))
ax1 = plt.subplot(141)
ax2 = plt.subplot(142)
ax3 = plt.subplot(143)
ax4 = plt.subplot(144)
    
for i in range(Nt):
    vel += acc * dt/2.0 # (1/2) kick
    pos += vel * dt # drift
    acc = getAcc( vel, pos ) # update accelerations
    vel += acc * dt/2.0 # (1/2) kick
    t += dt # update time
    pos = bound(pos)
    
    
    if float(i)/50. == round(float(i)/50.):
        hist, bins = histogram(pos)
        plt.sca(ax1)
        plt.cla()
        #plt.tricontourf(pos[:,0],pos[:,1], vel[:,0], 15)
        plt.scatter(pos[:,0],pos[:,1],s=10,color='blue')
        plt.xlim([0, 1]); plt.ylim([0, 1])
        plt.title('t=%0.5f'%t)
        
        plt.sca(ax2)
        plt.cla()
        counts,xbins,ybins=np.histogram2d(pos[:,0],pos[:,1],bins=30)
        plt.imshow(counts.T, interpolation='bilinear', origin='lower',\
        extent=[0,1,0,1],vmin=0,vmax=20)
        plt.xlim([0, 1]); plt.ylim([0, 1])
        plt.title('density map')
        
        plt.sca(ax3)
        plt.cla()
        plt.plot(bins,hist)
        plt.xlabel('y')
        plt.ylabel('particle number')
        plt.title('X-averaged density')
        plt.ylim([100, 300])
        
        plt.sca(ax4)
        plt.cla()
        plt.tricontourf(pos[:,0],pos[:,1], (acc[:,0]), 15,vmax=100,vmin=-100)
        plt.xlim([0, 1]); plt.ylim([0, 1])
        plt.title('$v_x$')
        
        #plt.savefig('fig%03d.png'%(float(i)/50))
        plt.show()
        
