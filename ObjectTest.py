import numpy as np
import random as r
import math as m
import matplotlib.pyplot as plt
import time as t

npartbase=20
hbase=0.01
gbase=0
boxbase=np.vstack([0,100])
plistposbase=[]
plistvelbase=[]
partbase={
    "spring":250,
    "radius":0.2
}

for n in range(npartbase):
    plistposbase.append([r.random()*100,r.random()*100])
    plistvelbase.append([r.random()*1,r.random()*1])

plistposarray=np.array(plistposbase)

def SimulationStep(plistpos=plistposarray,plistvel=plistvelbase,h=hbase,part=partbase,g=gbase,npart=npartbase):
    Forces=[]
    for i in range(npart):
        Forces.append(0)
        xnew=plistpos[i][0]+h*plistvel[i][0]+h**2*Forces[i]
        ynew=plistpos[i][1]+h*plistvel[i][1]+h**2*Forces[i]
        vxnew=(xnew-plistpos[i][0])/h
        vynew=(ynew-plistpos[i][1])/h
    plt.scatter(*plistpos.T)
    
    
    return(h)

for b in range(10):
    SimulationStep()
    t.sleep(1)
