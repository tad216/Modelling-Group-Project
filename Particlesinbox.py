# setup packages
import numpy as np
import matplotlib.pyplot as plt
import math

n = 50 # number of particles
box = np.array([[0,10*np.sqrt(n)],[0,10*np.sqrt(n)]]) # 2x2 array showing size of box

p = np.zeros((n,2)) # empty position array
xrand = np.random.rand(50,1) # random distribution of x-points from 0-1
yrand = np.random.rand(50,1) # y-points
xrand = np.round(xrand*box[0,1],5) # fit the positions to the box size and round them slightly
yrand = np.round(yrand*box[1,1],5) 
for m in range(50): # replace empty positions with randomised values
    p[m,0] = xrand[m,0]
    p[m,1] = yrand[m,0]
    
v = np.zeros((n,2)) # same thing for starting vectors
vini = 2.5 # suggested intitial speed
vxrand = np.random.rand(50,1) 
vyrand = np.random.rand(50,1) 
vxrand = np.round(2*np.subtract(vxrand,0.5)*vini,5) # randomise the vector sign and fit to bounds of vini
vyrand = np.round(2*np.subtract(vyrand,0.5)*vini,5) 
for m in range(50): # replace empty vector positions
    v[m,0] = vxrand[m,0]
    v[m,1] = vyrand[m,0] 


def d(p1,p2): # find distance between particles
    return(np.sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2))
def a(p1,p2): # angle between points 
    return(math.atan((p1[1]-p2[1])/(p1[0]-p2[0])))
    
def Simulationstep(p,v,h,r,k,box,g): # time-step function
    f = np.zeros((n,2)) # empty array of forces at step n
    for pm in range(n): # update forces array, checking at point pm
        for pnotm in range(pm+1,n): # check every other point from then onwards, not counting backwards
            dnow=d(p[pm],p[pnotm])
            if dnow<2*r: # if particles within range
                anow=a(p[pm],p[pnotm])
                f[pm,0] = k*(2*r-dnow)*math.cos(anow) # set new force, replacing zeros
                f[pm,1] = k*(2*r-dnow)*math.sin(anow)
            else:
                continue
            
        fleft = max(0, r+box[0,0]-p[pm,0]) # force at left wall is either positive or 0
        fright = max(0, r-box[0,1]+p[pm,0]) 
        fbot = max(0, r+box[1,0]-p[pm,1])
        ftop = max(0, r-box[1,1]+p[pm,1])
        fwall = np.array([[k*(fleft-fright),k*(fbot-ftop)]]) # vector of wall forces on point pm
        f[pm]=np.add(f[pm],fwall[0]) # add wall force on point pm to current force (probably 0)
        
        # gravity calculation goes here
    
    pnew = np.add(p, h*v, h**2*f) # use verlet formula with the 3 arrays position, velocity, force
    v = np.subtract(pnew, p)/h # velocity=displacement/time
    p = pnew
    return(p,v)


h = 0.01 # time step
r = 0.2 # radius of particle
k = 250 # spring constant
g = 0 # acceleration due to gravity
                
    