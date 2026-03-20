# Import Functions
import numpy as np
import random as r
import math as m
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time as t
import copy as copy

# Set up parameters
npartbase=4**5 #Number of particles
hbase=0.01 # Timestep
gbase=0 #Gravity (positive values increase gravity)
def initialiseVariables(npartbase=npartbase):  
    global boxbase,plistposbase,plistvelbase,partbase,l,u,densbase,colours,Temp,dp,vnew_0,wc,pc,pc_vals,wc_vals,iter,n_t,n_rec,temp_vals,t_1
    boxbase=10*(npartbase**0.5) #Box size
    #Dictionary for constants
    partbase={
    
        "spring":250,
        "radius":1
    }
    plistposbase=np.zeros((npartbase,2)) # Empty position
    plistvelbase=np.zeros((npartbase,2)) # Empty velocity
    
    # Create arrays of initial velocities and positions for the particles
    for n in range(npartbase):
        plistposbase[n,:]=([r.random()*boxbase,r.random()*boxbase]) #Initial positions
        plistvelbase[n,:]=([(r.random()-0.5)*2*2.5,(r.random()-0.5)*2*2.5]) # Initial Velocities (Originally 2.5)
    
    l = np.array([0,0]) #Coordinates of the bottom left corner of the box
    u = np.array([boxbase,boxbase]) #Coordinates of the top right corner of the box
    densbase = npartbase/((u[0]-l[0])*(u[1]-l[1])) #Density of particles
    
    # Give each particle a unique colour
    colours=[]
    for n in range(npartbase):
        colours.append([n/npartbase,0,1-n/npartbase,1])
    
    #Here we intitialize all the arrays that we will need for our functions
    Temp = np.zeros(npartbase) #Empty Temperature array
    temp_vals = []
    dp = np.zeros([npartbase]) #Change in distance per timestep
    vnew_0 = np.zeros([npartbase]) #RMS of velocities
    wc = np.zeros([npartbase,2]) #Confirmation of a wall collision
    pc = np.zeros([npartbase,npartbase]) #Confirmation of a particle collision
    pc_vals = np.zeros([npartbase]) #Number of particle collisions
    wc_vals = np.zeros([npartbase]) #Number of wall collisions
    iter=0 #Number of iterations
    n_t = 1000
    n_rec = int(n_t/2)
    t_1 = 0
    a = 2

initialiseVariables() 

def wall_coll(i, p,part,npart = npartbase, iter = iter): #Collisions with the wall function
    global wc
    global wc_vals
    #The force exerted on the particle is proposional to how far inside the wall the particle is
    f_left = max(0,part["radius"] + l[0]- p[i,0])
    f_right = max(0,part["radius"]+ p[i,0] - u[0])
    f_bot = max(0,part["radius"] + l[1] - p[i,1])
    f_up = max(0, part["radius"] + p[i,1] - u[1])
    #The spring constant affects how long the collision takes. the higher the spring constant the shorter the collision
    F_C = part["spring"]*np.array([f_left - f_right, f_bot - f_up])
    #This if statement makes collisions only recorded after a certain amount of iterations
    if iter>n_rec:
        if f_left !=0 or f_right !=0: #Left and Right wall collision tracker
                if wc[i,0] == 0:
                    wc[i,0] = 1
                    #A collision is only counted if there was not a collision on the last step,
                    #and only ends when there is no collision on the current step.
                    #This ensures that each collision is counted only once.
                    wc_vals[i] += 1
                    #print("side collision",wc_vals)
        else: 
            wc[i,0] = 0
        if f_up != 0 or f_bot !=0: #Top and Bottom wall collision tracker
            if wc[i,1] == 0:
                    wc[i,1] = 1
                    #This functions the same as the above if statement
                    wc_vals[i] += 1
                    #print("top/bottom collision",wc_vals)
        else: wc[i,1] = 0
    else: wc_vals[i] = 0
    return F_C



def particle_coll(i,j,p,part,npart,iter = iter):  #Collisions with other particles function
    global pc
    global pc_vals
    #This 2d array holds the x and y components of the force between the particles
    force = np.array([0,0])
    d = np.sqrt((p[i,0]-p[j,0])**2 +(p[i,1]-p[j,1])**2)
    alpha = np.arctan2(p[i,1]-p[j,1],p[i,0]-p[j,0])
    #We check to make sure that the particles are colliding here by checking the distance between them
    if 0 < d < 2*part["radius"]:
        #The force between particles is calculated the same as the wall, using the spring constant
        force = part["spring"]*(2*part["radius"] - d)*np.array([np.cos(alpha),np.sin(alpha)])
    #This if statement has the same purpose as the one in the wall collision function
    if iter>n_rec:
        if 0 < d < 2*part["radius"]:
                if pc[i,j] == 0:
                    pc[i,j] = 1
                    pc_vals[i] += 1
                    #print("particle collision", pc_vals)
        else: pc[i,j] == 0
    else: pc_vals[i] = 0
    return force

# define function to combine p and v
def combine(array1, array2):
    combined = np.concatenate((array1, array2), axis=1)
    return(combined)
def split(array):
    split1 = np.delete(array,np.s_[2,3],1)
    split2 = np.delete(array,np.s_[0,1],1)
    return(split1, split2)

#This is the function that gets called to perform each step
def SimulationStep(p=plistposbase,v=plistvelbase,h=hbase,part=partbase,g=gbase,npart=npartbase,dp = dp, vnew_0 = vnew_0,iter = iter):
    combined = combine(p,v)
    combined = combined[combined[:, 0].argsort()]
    #Some local arrays are initiated
    F = np.zeros((npart,2))
    Temp = np.zeros(npart)
    p,v=split(combined)
    #This part runs the collision functions, for each pair of particles
    for i in range(npart):
        for j in range(i+1, npart):
            if (np.abs(p[j,0]-p[i,0])) > 2*part["radius"]: # stop for loop if particles are a minimum of 2r apart
                break
            else:
                force = particle_coll(i,j,p,part,npart,iter = iter)     
                F[i,:] += force
                F[j,:] += -force
        F_C = wall_coll(i,p,part,npart,iter = iter)
        F[i,:] += F_C
            # Includes Gravity
        F[i,1] += -g
    # The positions of the particles are changed using the Verlet updating formula 
    pnew = p + h*v + (h**2)*F
    vnew =(pnew - p)/h
    # Temperature
    mean_temp_old = np.mean(Temp)
    for i in range(npart):
        Temp[i] = 0.5*np.sum(vnew[i,:]**2)
    mean_temp = np.mean(Temp)
    for i in range(npart):
        Temp[i] = 0.5*np.sum(vnew[i,:]**2)
    mean_temp = np.mean(Temp)
    temp_vals.append(mean_temp)
    temp_av = np.mean(temp_vals)
    for m in range(npartbase):
        vnew_0[m] = np.sqrt((vnew[m,0])**2 +(vnew[m,1])**2)
    if iter>n_rec:
        dp += vnew_0*h
    else: 
         dp += np.zeros([npartbase])
    #print("action", action)
    #print("vnew",vnew)
    
    return(pnew,vnew, Temp, mean_temp, mean_temp_old, dp,vnew_0)


begin = t.time()

# Create plot for animating the particles
def animation():
    global t_1
    plt.ion()
    fig=plt.figure()
    ax=fig.add_subplot(111)
    plpts = ax.scatter(*plistposbase.T,c=colours)
    plt.xlim(0,boxbase)
    plt.ylim(0,boxbase)
    xnext=plistposbase
    vnext=plistvelbase
    input("press enter")
    for iter in range(n_t):
        #Saves the new x and v for the next step, and the temps and path for data collection
        xnext,vnext, Temp, mean_temp,mean_temp_old,dp,vnew_0=SimulationStep(xnext,vnext,iter = iter)
        if np.abs(mean_temp - mean_temp_old) < 0.1:
            end = t.time()
            t_1 = end - begin
            t_2 = t_1 +5*np.sqrt(npartbase)/a
        #This updates the graph
        plpts.set_offsets(xnext)
        plt.pause(0.01)
        plt.show()
    print(t_1)


def sim(plistposbase = plistposbase, plistvelbase = plistvelbase):
    global t_1
    global vnew_vals
    for iter in range(n_t):
        plistposbase, plistvelbase, Temp, mean_temp, mean_temp_old,dp,vnew_0,=SimulationStep(plistposbase,plistvelbase,iter = iter)
        print(iter)
        
        
        
#print(xnext)
# Imput to space out code

def func():
    print("wc", wc_vals,"pc", pc_vals)
    mean_path = dp/(1+wc_vals+pc_vals)
    print("path",mean_path)
    check = np.sum(mean_path)/npartbase
    dp_c = np.mean(dp)
    print("Average Distance Travelled",dp_c)
    print("mean path per collision", check)



sim()


#not the complete thing, update if you can
"""pmax = 4
N_vals = []
stand_dev = []
for j in range(pmax):
    N = 4**j
    initialiseVariables(N)
    Temp_vals = []
    for e in range(n_t):
        plistposbase,plistvelbase,Temp, mean_temp,dp = SimulationStep(plistposbase,plistvelbase, h=hbase, part = partbase, g=gbase , npart = N)
        if (e>n_rec):
            Temp_vals.append(mean_temp)
    N_vals.append(N)
    sd = np.std(Temp_vals[-n_rec-1:])
    stand_dev.append(sd)
    print("standard deviation", sd)
    #print(Temp_vals)
plt.loglog(N_vals,stand_dev,"*")

#print("Standard Dev", stand_dev)
print("v", vnew_0)
input("press enter")"""
plt.hist(vnew_0, bins = 32,edgecolor = "black") #x-axis = speed, y-axis = frequency 
plt.xlabel("Speed of Particles")
plt.ylabel("Frequency")
plt.show()
plt.savefig()
#print("Standard Dev", stand_dev)"""
