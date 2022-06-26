import numpy as np
import matplotlib.pyplot as plt
import sys



def initial_conditions(n, width, length):
    a=np.random.random((3,n))
    for i in range(np.ma.size(a, 1)):
        a[0,i]= int(np.ceil(2*a[0,i])) 
        if a[0,i]==1.0:
            a[0,i]=1.0
        elif a[0,i]==2.0:
            a[0,i]=-1.0    
        a[1,i]=(width-2.0)*a[1, i]+1.0
        a[2,i]=np.random.normal(1.3, 0.26) 
   
    old=np.zeros((4,n))
    for i in range(n):
        if a[0,i]==1.0:
            old[0,i]=0.0 #Initial x position
        elif a[0,i]==-1.0:
            old[0,i]=length #Initial x position   
        old[1, i]=a[1, i] #Initial y position
        old[2, i]=a[2, i]*a[0, i] #Initial x  velocity vx
        old[3, i]=0.0 #Initial y velocity vy
           
    return a, old   



def position_evolver(a, old, i, delta_t, tau, width, max_v_coef, pos_int_wall, widthwall, sigma, v0, r, u0, c, phi):
    n=np.ma.size(a,1)
    new=np.zeros(4)
    rabx=np.zeros(n)
    raby=np.zeros(n)
    fy=np.zeros(n)
    fx=np.zeros(n)
    for j in range(n):
        rabx[j]=old[0, i]-old[0, j]
        raby[j]=old[1, i]-old[1, j]
    #Goal force
    new[2]=delta_t*(1.0/tau*(a[2, i]*a[0, i]-old[2, i])) #Term 1
    new[3]=delta_t*(1.0/tau*(-old[3, i]))
    deltat2=delta_t
    #Pedestrian-pedestrian repulsion force
    for j in range(n):
        if j!=i:
      
            deltax=0.001
            rabmod=np.sqrt(rabx[j]**2.0+raby[j]**2.0)
            rabmodx=np.sqrt((rabx[j]+deltax)**2.0+raby[j]**2.0)
            rabmody=np.sqrt(rabx[j]**2.0+(raby[j]+deltax)**2.0)

            
            theta=np.arctan2(raby[j],rabx[j])
            thetax=np.arctan2(raby[j],rabx[j]+deltax)
            thetay=np.arctan2(raby[j]+deltax,rabx[j])
            vb=np.sqrt(old[2,j]**2.0+old[3, j]**2.0)
            root=np.sqrt(rabmod**2.0-2.0*vb*deltat2*rabmod*np.cos(theta)+vb**2.0*deltat2**2.0)
            rootx=np.sqrt(rabmodx**2.0-2.0*vb*deltat2*rabmodx*np.cos(thetax)+vb**2.0*deltat2**2.0)
            rooty=np.sqrt(rabmody**2.0-2.0*vb*deltat2*rabmody*np.cos(thetay)+vb**2.0*deltat2**2.0)
            b=np.sqrt(rabmod**2.0+2.0*rabmod*root+root**2.0-vb**2.0*deltat2**2.0)/2.0
            bx=np.sqrt(rabmodx**2.0+2.0*rabmodx*rootx+rootx**2.0-vb**2.0*deltat2**2.0)/2.0
            by=np.sqrt(rabmody**2.0+2.0*rabmody*rooty+rooty**2.0-vb**2.0*deltat2**2.0)/2.0
            exp=np.e**(-b/sigma)
            expx=np.e**(-bx/sigma)
            expy=np.e**(-by/sigma)
            fx[j]=-v0*(expx-exp)/deltax
            fy[j]=-v0*(expy-exp)/deltax
        
    
            if (-a[0,i]*fx[j]>=np.sqrt(fx[j]**2.0+fy[j]**2.0)*np.cos(phi)):
                new[2]+=delta_t*fx[j]
                new[3]+=delta_t*fy[j]
            else:
                new[2]+=delta_t*fx[j]*c
                new[3]+=delta_t*fy[j]*c

    distInf=np.abs(old[1, i])
    distSup=np.abs(width-old[1, i])    
 
    wallsInfy=u0/r*np.e**((-distInf)/r)
    wallsInfx=0.0
  
    
    new[2]+=delta_t*wallsInfx    
    new[3]+=delta_t*(wallsInfy)
    
    #Pedestrian-wall repulsion force
    wallsSupy=-u0/r*np.e**((-distSup)/r) 
    wallsSupx=0.0

    new[2]+=delta_t*wallsSupx          
    new[3]+=delta_t*wallsSupy  
    
    new[2]+=old[2, i]
    new[3]+=old[3, i]
    
    if np.sqrt(new[2]**2.0+new[3]**2.0)>max_v_coef*a[2, i]:
        new[2]=new[2]*max_v_coef*a[2, i]/np.sqrt(new[2]**2.0+new[3]**2.0)
        new[3]=new[3]*max_v_coef*a[2, i]/np.sqrt(new[2]**2.0+new[3]**2.0)
   

    new[0]=delta_t*new[2]+old[0,i]
    new[1]=delta_t*new[3]+old[1, i]
    
    
    #El problema de no incluir esto es que si no se pone la restricción sobre el punto inicial es probable
    #que se abandone por los lados inferior/superior los bordes
    if new[1]<=0.0 or new[1]>=width:
        new[0]=old[0,i]
        new[1]=old[1, i]

    

            
    
            
    return  new  







t_delta=0.1 #Have to mantain low deltat. From a physical point of view could be understood as the reaction time of a person
#ie: the time a person needs to analyse a situation and take a decision os their velocity and direction. From a computational
# point of view because the method used to solve is finite differences it cannot be very big, for t=0.1s works fine
t_max=52.0
t_step=int(t_max/t_delta)
n=150
width=10.0
length=50.0
tau=0.5
r=0.2
u0=10.0
max_v_coef=1.3
pos_int_wall=length/2.0
widthwall=2.0
sigma=0.3
v0=2.1
c=0.5
phi=100.0*2.0*np.pi/360



init, old=initial_conditions(n, width, length)

new=np.zeros((4, n))

for l in range(t_step):
    for i in range(n):
        new[:,i]=position_evolver(init, old, i, t_delta, tau, width, max_v_coef, pos_int_wall, widthwall, sigma, v0, r, u0, c, phi)
        
    if ((l*t_delta)%1.0)==0.0:
        azules=[]
        rojos=[]
        for k in range(n):
            if init[0, k]==1.0:
                azules.append([new[0,k], new[1, k],10*np.sqrt(new[2, k]**2.0+new[3,k]**2.0)])
            elif init[0,k]==-1.0:
                rojos.append([new[0,k], new[1, k], 10*np.sqrt(new[2, k]**2.0+new[3,k]**2.0)])
        nprojos=np.array(rojos)
        npazules=np.array(azules)

        try:
            plt.scatter(npazules[:,0], npazules[:,1],s=npazules[:,2],  c='b')
        except:
            pass
        try:
            plt.scatter(nprojos[:,0], nprojos[:,1],s=nprojos[:,2], c='r')
        except:
            pass
        plt.xlim(0,length)
        plt.ylim(0,width)
        plt.title(str(l*t_delta))
        plt.savefig('gif'+str(l*t_delta)+'.png')
        plt.show()
        print(l*t_delta)
        
    for i in range(4):
        for j in range(n):
            old[i,j]=new[i,j]    
     