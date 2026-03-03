import numpy as np
import random as r
import math as m
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time as t

npartbase=40
hbase=0.01
gbase=0
boxbase=10*(npartbase**0.5)
plistposbase=[]
plistvelbase=[]
partbase={
    "spring":250,
    "radius":0.4
}

plistvelbase=np.zeros((npartbase,2))
plistposbase=np.zeros((npartbase,2))

for n in range(npartbase):
    plistposbase[n,:]=([r.random()*boxbase,r.random()*boxbase])
    if (plistposbase[n,0]<partbase["radius"]):
        plistposbase[n,0]=partbase["radius"]+plistposbase[n,0]
    elif(plistposbase[n,0]>boxbase-partbase["radius"]):
        plistposbase[n,0]=plistposbase[n,0]-partbase["radius"]
    if (plistposbase[n,1]<partbase["radius"]):
        plistposbase[n,1]=partbase["radius"]+plistposbase[n,1]
    elif (plistposbase[n,1]>boxbase-partbase["radius"]):
        plistposbase[n,1]=plistposbase[n,1]-partbase["radius"]
    plistvelbase[n,:]=([(r.random()-0.5)*40,(r.random()-0.5)*40])

print(plistposbase,plistvelbase)

def SimulationStep(p=plistposbase,v=plistvelbase,h=hbase,part=partbase,g=gbase,npart=npartbase):
    """
    Forces=[]
    plistposnew=np.zeros((npart,2))
    plistvelnew=np.zeros((npart,2))
    for i in range(npart):
        Forces.append(0)
        xnew=plistpos[i,0]+h*plistvel[i,0]+h**2*Forces[i]
        ynew=plistpos[i,1]+h*plistvel[i,1]+h**2*Forces[i] 
        vxnew=(xnew-plistpos[i,0])/h
        vynew=(ynew-plistpos[i,1])/h
        if (xnew<partbase["radius"] or xnew>boxbase-partbase["radius"]):
            vxnew=-vxnew
        if (ynew<partbase["radius"] or ynew>boxbase-partbase["radius"]):
            vynew=-vynew
        for j in range(npart):
            if ((xnew-plistpos[j,0])**2+(ynew-plistpos[j,1])**2<=(partbase["radius"]*2)**2):
                velstore=[vxnew,vynew]
                vxnew=plistvel[j,0]
                vynew=plistvel[j,1]
                plistvel[j]=velstore
        
        plistposnew[i,:]=([xnew,ynew])
        plistvelnew[i,:]=([vxnew,vynew])
    """
    F = np.zeros((npart,2))
    for i in range(npart):
        for j in range(i+1, npart):
            if i!=j:
                d = np.sqrt((p[i,0]-p[j,0])**2 +(p[i,1]-p[j,1])**2)
                alpha = np.arctan2(p[i,1]-p[j,1],p[i,0]-p[j,0])
            else:
                 d = 0
            if 0 < d < 2*part["radius"]:
                 force = part["spring"]*(2*part["radius"] - d)*np.array([np.cos(alpha),np.sin(alpha)])
                 F[i,:] += force
                 F[j,:] += -force
        f_left = max(0,part["radius"] - p[i,0])
        f_right = max(0,part["radius"]+ p[i,0] - 10*np.sqrt(npart))
        f_bot = max(0,part["radius"]- p[i,1])
        f_up = max(0, part["radius"] + p[i,1] - 10*np.sqrt(npart))
        F_C = part["spring"]*np.array([f_left - f_right, f_bot - f_up])
        F[i,:] += F_C
        #print(F)
    pnew = p + h*v + (h**2)*F
    vnew =(pnew - p)/h
    
    return(pnew,vnew)

plt.ion()
fig=plt.figure()
ax=fig.add_subplot(111)
plpts = ax.scatter(*plistposbase.T)
plt.xlim(0,boxbase)
plt.ylim(0,boxbase)
xnext=plistposbase
vnext=plistvelbase

for b in range(500):
    xnext,vnext=SimulationStep(xnext,vnext)
    #print(xnext[0][0])
    plpts.set_offsets(xnext)
    plt.pause(0.01)
    plt.show()

input("Press enter")
