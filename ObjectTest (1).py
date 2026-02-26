import numpy as np
import random as r
import math as m
import matplotlib.pyplot as plt
import matplotlib.animation as animation
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
    if (plistposbase[n][0]<partbase["radius"]):
        plistposbase[n][0]=partbase["radius"]+plistposbase[n][0]
    elif(plistposbase[n][0]>100-partbase["radius"]):
        plistposbase[n][0]=plistposbase[n][0]-partbase["radius"]
    if (plistposbase[n][1]<partbase["radius"]):
        plistposbase[n][1]=partbase["radius"]+plistposbase[n][1]
    elif (plistposbase[n][1]>100-partbase["radius"]):
        plistposbase[n][1]=plistposbase[n][1]-partbase["radius"]
    plistvelbase.append([(r.random()-0.5)*200,(r.random()-0.5)*200])

plistposarray=np.array(plistposbase)

def SimulationStep(plistpos=plistposarray,plistvel=plistvelbase,h=hbase,part=partbase,g=gbase,npart=npartbase):
    Forces=[]
    plistposnew=[]
    plistvelnew=[]
    for i in range(npart):
        Forces.append(0)
        xnew=plistpos[i][0]+h*plistvel[i][0]+h**2*Forces[i]
        ynew=plistpos[i][1]+h*plistvel[i][1]+h**2*Forces[i] 
        vxnew=(xnew-plistpos[i][0])/h
        vynew=(ynew-plistpos[i][1])/h
        if (xnew<partbase["radius"] or xnew>100-partbase["radius"]):
            vxnew=-vxnew
        if (ynew<partbase["radius"] or ynew>100-partbase["radius"]):
            vynew=-vynew
        for j in range(npart):
            if ((xnew-plistpos[j][0])**2+(ynew-plistpos[j][1])**2<=(partbase["radius"]*2)**2):
                velstore=[vxnew,vynew]
                vxnew=plistvel[j][0]
                vynew=plistvel[j][1]
                plistvel[j]=velstore
        
        plistposnew.append([xnew,ynew])
        plistvelnew.append([vxnew,vynew])
    
    return(plistposnew,plistvelnew)

plt.ion()
fig=plt.figure()
ax=fig.add_subplot(111)
plpts = ax.scatter(*plistposarray.T)
plt.xlim(0,100)
plt.ylim(0,100)
xnext=plistposbase
vnext=plistvelbase

for b in range(1000):
    xnext,vnext=SimulationStep(xnext,vnext)
    print(xnext[0][0])
    plpts.set_offsets(xnext)
    plt.pause(0.01)
    plt.show()

input("Press enter")

