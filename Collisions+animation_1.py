import numpy as np
import matplotlib.pyplot as plt
rand = np.random.rand
from math import *

p = 3
N = 4**p # number of particles
K = 250 # spring constant
r = 0.2 # radius of particles
h = 0.01 # time step size

l_1 = 0
l_2 = 0
l = np.array([l_1,l_2]) # lower-left corner of the box
u_1 = 10*np.sqrt(N)
u_2 = 10*np.sqrt(N)
u = np.array([u_1,u_2]) # upper-right corner of the box
g = 0 # gravity
tini = 0
tend = 20
t = 0
xini = rand(2,N)
vini = 5

x = np.vstack([l[0]+rand(1,N)*(u[0]-l[0]), 
               l[1]+rand(1,N)*(u[1]-l[1])])

v = 2*(rand(2,N)-0.5)*vini

part = list("spring" = K, "radius" = r)
box = list(l_1,l_2,u_1,u_2)



for  m in range(100):
    F = np.zeros((2,N))
    for i in range(N):
        for j in range(i+1, N):
            if i!=j:
                d = np.sqrt((x[0,i]-x[0,j])**2 +(x[1,i]-x[1,j])**2)
                alpha = np.arctan2(x[1,i]-x[1,j],x[0,i]-x[0,j])
            else:
                 d = 0
            if 0 < d < 2*r:
                 force = K*(2*r - d)*np.array([np.cos(alpha),np.sin(alpha)])
                 F[:,i] += force
                 F[:,j] += -force
        f_left = max(0,r+l_1 - x[0,i])
        f_right = max(0,r+ x[0,i] - 10*np.sqrt(N))
        f_bot = max(0,r+ l_2 - x[1,i])
        f_up = max(0, r + x[1,i] - 10*np.sqrt(N))
        F_C = K*np.array([f_left - f_right, f_bot - f_up])
        F[:,i] += F_C
    x_new = x + h*v + (h**2)*F
    v_new =(x_new - x)/h
    v = v_new
    x = x_new
    t = t+h
    plt.plot([l[0],u[0], u[0],l[0],l[0]],
        [l[1],l[1],u[1],u[1],l[1]])
    plt.plot(x[0,:],x[1,:],".")
    plt.pause(0.01)
    plt.clf()

plt.show

plt.show
