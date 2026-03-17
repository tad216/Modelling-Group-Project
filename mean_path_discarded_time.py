# Import Functions
import numpy as np
import random as r
import math as m
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time as t

# Set up parameters
npartbase=32 #Number of particles
hbase=0.01 # Timestep
gbase=0 #Gravity (positive values increase gravity)
boxbase=10*(npartbase**0.5) #Box size
plistposbase=[] # Positions
plistvelbase=[] # Velocity
#Dictionary for constants
partbase={

    "spring":250,
    "radius":0.2
}
plistposbase=np.zeros((npartbase,2)) # Empty position
plistvelbase=np.zeros((npartbase,2)) # Empty velocity

# Create arrays of initial velocities and positions for the particles
for n in range(npartbase):
    plistposbase[n,:]=([r.random()*boxbase,r.random()*boxbase]) #Initial positions
    plistvelbase[n,:]=([(r.random()-0.5)*2.5,(r.random()-0.5)*2.5]) # Initial Velocities (Originally 2.5)

l = np.array([0,0]) #Coordinates of the bottom left corner of the box
u = np.array([boxbase,boxbase]) #Coordinates of the top right corner of the box

densbase = npartbase/((u[0,0]-l[0,0])*(u[0,1]-l[0,1])) #Density of particles

# Give each particle a unique colour
colours=[]
for n in range(npartbase):
    colours.append([r.random(),r.random(),r.random(),1])


Temp = np.zeros(npartbase) #Empty Temperature array
dp = np.zeros([npartbase]) #Change in distance per timestep
vnew_0 = np.zeros([npartbase]) #RMS of velocities

wc = np.zeros([npartbase,2]) #Confirmation of a wall collision
pc = np.zeros([npartbase,npartbase]) #Confirmation of a particle collision
pc_vals = np.zeros([npartbase]) #Number of particle collisions
wc_vals = np.zeros([npartbase]) #Number of wall collisions
iter=0 #Number of iterations 

def wall_coll(i, p,part,npart = npartbase, iter = iter): #Collisions with the wall function
    global wc
    global wc_vals
    f_left = max(0,part["radius"] + l[0,0]- p[i,0])
    f_right = max(0,part["radius"]+ p[i,0] - u[0,0])
    f_bot = max(0,part["radius"] + l[0,1] - p[i,1])
    f_up = max(0, part["radius"] + p[i,1] - u[0,1])
    F_C = part["spring"]*np.array([f_left - f_right, f_bot - f_up])
    if iter>n_rec:
        if f_left !=0 or f_right !=0: #Left and Right wall collision tracker
                if wc[i,0] == 0:
                    wc[i,0] = 1
                    wc_vals[i] += 1
                    print("side collision",wc_vals)
        else: 
            wc[i,0] = 0
        if f_up != 0 or f_bot !=0: #Top and Bottom wall collision tracker
            if wc[i,1] == 0:
                    wc[i,1] = 1
                    wc_vals[i] += 1
                    print("top/bottom collision",wc_vals)
        else: wc[i,1] = 0
    else: wc_vals[i] = 0
    return F_C



def particle_coll(i,j,p,part,npart,iter = iter):  #Collisions with other particles function
    global pc
    global pc_vals
    force = np.array([0,0])
    if i!=j:
        d = np.sqrt((p[i,0]-p[j,0])**2 +(p[i,1]-p[j,1])**2)
        alpha = np.arctan2(p[i,1]-p[j,1],p[i,0]-p[j,0])
    else:
        d = 0
    if 0 < d < 2*part["radius"]:
        force = part["spring"]*(2*part["radius"] - d)*np.array([np.cos(alpha),np.sin(alpha)])
    if iter>n_rec:
        if 0 < d < 2*part["radius"]:
                if pc[i,j] == 0:
                    pc[i,j] = 1
                    pc_vals[i] += 1
                    print("particle collision", pc_vals)
        else: pc[i,j] == 0
    else: pc_vals[i] = 0
    return force


#print(plistposbase,plistvelbase)

def SimulationStep(p=plistposbase,v=plistvelbase,h=hbase,part=partbase,g=gbase,npart=npartbase,dp = dp, vnew_0 = vnew_0,iter = iter):
    F = np.zeros((npart,2))
    Temp = np.zeros(npart)
    for i in range(npart):
        for j in range(i+1, npart):
            force = particle_coll(i,j,p,part,npart,iter = iter)     
            F[i,:] += force
            F[j,:] += -force
        F_C = wall_coll(i,p,part,npart,iter = iter)
        F[i,:] += F_C
        F[i,1] += -g
    #print("Forces",F) 
    # verlet updating formula 
    pnew = p + h*v + (h**2)*F
    vnew =(pnew - p)/h
    # Temperature
    for i in range(npart):
        Temp[i] = 0.5*np.sum(vnew[i,:]**2)
    mean_temp = np.mean(Temp)
    if iter>n_rec:
        vnew_0 = np.sqrt((vnew[:,0])**2 +(vnew[:,1])**2)
        dp += vnew_0*h
    else: 
         vnew_0 = np.zeros([npartbase])
         dp +=vnew_0*h

    #print("dp",dp)
    print(iter)
    
    return(pnew,vnew, Temp, mean_temp, dp)


#SimulationStep()
# Create box for particles
plt.ion()
fig=plt.figure()
ax=fig.add_subplot(111)
plpts = ax.scatter(*plistposbase.T,c=colours)
plt.xlim(0,boxbase)
plt.ylim(0,boxbase)
xnext=plistposbase
vnext=plistvelbase
input("press enter")

n_t = 1000
n_rec = int(n_t/2)

# Animation Loop
for iter in range(n_t):
    xnext,vnext, Temp, mean_temp,dp=SimulationStep(xnext,vnext,iter = iter)
    #print(xnext[0][0])
    plpts.set_offsets(xnext)
    plt.pause(0.01)
    plt.show()
#print(xnext)
# Imput to space out code
input("Press enter")
print("wc", wc_vals,"pc", pc_vals)
mean_path = dp/(1+wc_vals+pc_vals)
print("path",mean_path)
check = np.sum(mean_path)/npartbase
dp_c = np.mean(dp)
print("Average Distance Travelled",dp_c)
print("mean path per collision", check)



#not the complete thing, update if you can
"""pmax = 6
N_vals = []
stand_dev = []
for j in range(pmax):
    N = 4**j
    vini = 2.5
    p = np.random.rand(N,2)
    v = 2*(np.random.rand(N,2)-0.5)*vini
    Temp_vals = []
    for e in range(n_t):
        p,v,Temp, mean_temp = SimulationStep(p = p, v =v, h=hbase, part = partbase, g=gbase , npart = N)
        Temp_vals.append(mean_temp)
    N_vals.append(N)
    sd = np.std(Temp_vals[-nrec-1:])
    stand_dev.append(sd)
    print("standard deviation", sd)
plt.loglog(N_vals,stand_dev,"*")

plt.show()

#print("Standard Dev", stand_dev)"""

print(densbase)
