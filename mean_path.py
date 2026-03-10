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
gbase=0 #
boxbase=10*(npartbase**0.5) #Boc size
plistposbase=[] # Positions
plistvelbase=[] # Velocity
#Dictionary for constants
partbase={

    "spring":250,
    "radius":1
}
plistposbase=np.zeros((npartbase,2)) # Empty position
plistvelbase=np.zeros((npartbase,2)) # Empty velocity
#Initial seed for testing code, 30 particles
"""
plistposbase=np.array([[32.77212199,28.37445098],
 [55.68867184,57.80481431],
 [52.84617268, 0.69400532],
 [50.82306856,59.21061824],
 [48.33877728, 5.95615546],
 [17.87698546,19.34803884],
 [ 6.83361335,55.25560096],
 [59.2023912 ,30.06644096],
 [36.01071364, 8.83980776],
 [34.42018804,45.18076791],
 [32.49656893,53.25594834],
 [49.7787714 ,25.83993974],
 [41.99174077,50.50557406],
 [55.63168024,22.32849789],
 [60.53838529,41.75819929],
 [46.6166008 ,33.55466903],
 [ 3.61000563,11.48461932],
 [17.65202004,16.80109725],
 [32.7185959 ,55.00030975],
 [17.32739128,49.48855923],
 [14.48216066,28.30471808],
 [ 0.42611602,28.14578894],
 [ 1.4944023 ,16.89600633],
 [28.50693849,27.53533577],
 [11.08823051,42.75261966],
 [45.16764202,21.35587934],
 [36.40537109,52.47699375],
 [31.17347544,27.56239692],
 [26.32771204,56.7481867 ],
 [16.00682906,62.56045091],
 [34.00746153,52.11333148],
 [ 9.7556969 ,20.70113218],
 [52.74736132, 2.6238682 ],
 [16.95141627,57.94668193],
 [25.1315872 ,35.43775712],
 [22.56269968,20.21882923],
 [61.9986072 ,23.37061282],
 [22.02443182,52.73742027],
 [61.52432831,14.30640276],
 [25.44575339,52.59073463]])
plistvelbase=np.array([[  6.91136697 ,11.12753943],
 [  9.84513795,-19.98124442],
 [ 16.01290807,   7.05480928],
 [  6.41783382, -7.79648871],
 [-17.01713444,  8.38253162],
 [ 19.57478054,  5.45015884],
 [ 10.55375394,-16.27086162],
 [-19.10522313,-17.90857429],
 [ -4.99441997, 15.7191499 ],
 [  5.63352979, 18.78125916],
 [ -7.86388004, 17.72610036],
 [ 19.95473693, -4.68249747],
 [ -6.50715986,-15.63394078],
 [ -4.04691294, 12.71092484],
 [ 13.13940514,  7.27873874],
 [-14.634454  ,-17.33511977],
 [ 15.47569635, -8.50532607],
 [ 14.49435857, -6.37929631],
 [  3.32208084,  8.48380402],
 [-18.3537866 ,  5.05094271],
 [ 16.4777879 ,-12.47035654],
 [-12.10481305,  6.42434027],
 [-16.61387156, -1.29733828],
 [ 14.38940133, 18.44601096],
 [-17.4570088 , -0.63369945],
 [ -7.8819144 , 10.69001836],
 [-17.49639382,-14.22333669],
 [  0.90050025,-11.36491525],
 [ -9.30885252, -4.99653324],
 [ 16.14858115,-17.2949247 ],
 [ 10.45877474, 11.42636387],
 [ -7.38116212, 15.92528399],
 [ -4.83558389,  7.50938653],
 [-14.04339256, -0.66882001],
 [  3.52055928, 18.41312515],
 [-19.38567694, 12.32947944],
 [-16.4705948 ,-11.6744364 ],
 [ 15.7874379 ,-19.65611136],
 [ -7.10448433,  2.00421728],
 [-11.82332228,-18.32302596]])
"""
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
    plistvelbase[n,:]=([(r.random()-0.5)*2.5,(r.random()-0.5)*2.5]) # Initial Velocities (Originally 2.5)

# Give each particle a unique colour
colours=[]
for n in range(npartbase):
    colours.append([r.random(),r.random(),r.random(),1])


Temp = np.zeros(npartbase) #Zero array of size N
dp = np.zeros([npartbase])
vnew_0 = np.zeros([npartbase])

wc = np.zeros([npartbase,2])
pc = np.zeros([npartbase,npartbase])
pc_vals = np.zeros([npartbase])
wc_vals = np.zeros([npartbase])


def wall_coll(i, p,part,npart): #Collisions with the wall function
    global wc
    global wc_vals
    f_left = max(0,part["radius"] - p[i,0])
    f_right = max(0,part["radius"]+ p[i,0] - 10*np.sqrt(npart))
    f_bot = max(0,part["radius"]- p[i,1])
    f_up = max(0, part["radius"] + p[i,1] - 10*np.sqrt(npart))
    F_C = part["spring"]*np.array([f_left - f_right, f_bot - f_up])
    if f_left !=0 or f_right !=0:
                if wc[i,0] == 0:
                    wc[i,0] = 1
                    wc_vals[i] += 1
                    print("side",wc_vals)
    else: 
        wc[i,0] = 0
    if f_up != 0 or f_bot !=0:
                if wc[i,1] == 0:
                    wc[i,1] = 1
                    wc_vals[i] += 1
                    print("top/bottom",wc_vals)
    else: wc[i,1] = 0
    
    return F_C



def particle_coll(i,j,p,part,npart):  #Collisions with other particles function
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
    if 0 < d < 2*part["radius"]:
                if pc[i,j] == 0:
                    pc[i,j] = 1
                    pc_vals[i] += 1
    else: pc[i,j] == 0
    return force


#print(plistposbase,plistvelbase)

def SimulationStep(p=plistposbase,v=plistvelbase,h=hbase,part=partbase,g=gbase,npart=npartbase,dp = dp, vnew_0 = vnew_0):
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
    Temp = np.zeros(npart)
    for i in range(npart):
        for j in range(i+1, npart):
            force = particle_coll(i,j,p,part,npart)     
            F[i,:] += force
            F[j,:] += -force
        F_C = wall_coll(i,p,part,npart)
        F[i,:] += F_C
    #print("Forces",F) 
    # verlet updating formula 
    pnew = p + h*v + (h**2)*F
    vnew =(pnew - p)/h
    # Temperature
    for i in range(npart):
        Temp[i] = 0.5*np.sum(vnew[i,:]**2)
    mean_temp = (1/npart)*np.sum(Temp)
    #print("Temp:", Temp, "Mean Temperature:", mean_temp)
    vnew_0 = np.sqrt((vnew[:,0])**2 +(vnew[:,1])**2)
    dp += vnew_0*h
    #print("dp",dp)
    
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

# Animation Loop
for b in range(200):
    xnext,vnext, Temp, mean_temp,dp=SimulationStep(xnext,vnext)
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
print("dp_c",dp_c)
print(check)
#not the complete thing, update if you can
"""pmax = 6
n_t = 100
nrec = int(n_t/2)
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

plt.show()"""

#print("Standard Dev", stand_dev)
