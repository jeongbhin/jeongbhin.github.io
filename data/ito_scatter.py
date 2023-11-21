# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 14:39:01 2023

@author: jbseo
"""

import numpy as np
import matplotlib.pyplot as plt


vx = 0
vy= 0
divV = 0.2


kappa_para = 1
kappa_perp = 0.00001*kappa_para

Bx = 1
By = 3
B = np.sqrt(Bx**2+By**2)

delkappa = 0.1
dt = 0.1
t = 0

X = 0; Y = 0; p = 10

alpha1 = np.array([np.sqrt(2*kappa_perp),0]).T
alpha2 = np.array([0,np.sqrt(2*kappa_perp)]).T
alpha3 = np.sqrt(2*(kappa_para-kappa_perp))*np.array([Bx/B,By/B]).T

X_sim =[] ; Y_sim = [] ; p_sim = []
for j in range(100):
    dW1, dW2, dW3 = np.sqrt(dt)*np.random.randn(3) 
    second_x = (alpha1[0]*dW1+alpha2[0]*dW2+alpha3[0]*dW3)
    second_y = (alpha1[1]*dW1+alpha2[1]*dW2+alpha3[1]*dW3)
    X += (delkappa+vx)*dt + second_x
    Y += (delkappa+vy)*dt + second_y
    p += -p/3*divV*dt
    X_sim.append(X) ; Y_sim.append(Y) ; p_sim.append(p)
    
    t += dt
    
x = np.linspace(min(X_sim),max(X_sim),100)    
plt.plot(x,By/Bx*x,color='r')
plt.plot(X_sim,Y_sim,color='k',alpha=0.5)    
plt.scatter(X_sim,Y_sim,c=p_sim)   
plt.colorbar()

