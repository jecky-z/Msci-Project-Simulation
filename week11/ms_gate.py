# -*- coding: utf-8 -*-
"""
Created on Sat Nov  6 19:45:33 2021

@author: jecky
"""

#%%
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import math
from qutip import *

options=qutip.solver.Options(atol=1e-07,rtol=1e-07,store_states=True)
params = {
   'axes.labelsize': 25,
   'font.size': 25,
   'legend.fontsize': 18,
   'xtick.labelsize': 25,
   'ytick.labelsize': 25,
   'figure.figsize': [20,10]
   } 
plt.rcParams.update(params)
###################################
name=('trap_frequency(2*pi), omega0(nu), eta, detuning(nu), size, initial state')
#############################
s_plus=(sigmax()+1j*sigmay())*0.5
s_minus=(sigmax()-1j*sigmay())*0.5

def H_sideband2(tlist,nu,omega0,eta,detuning,phase,size,order=5,form=0):
    if form ==0:# Array format
        c=0.5*hbar*omega0
        H=[]
        H.append([c*tensor(s_plus,identity(2),identity(size)),np.exp(1j*((detuning*tlist)+phase))])
        H.append([c*tensor(identity(2),s_plus,identity(size)),np.exp(1j*((detuning*tlist)+phase))])
        H.append([c*tensor(s_minus,identity(2),identity(size)),np.exp(-1j*((detuning*tlist)+phase))])
        H.append([c*tensor(identity(2),s_minus,identity(size)),np.exp(-1j*((detuning*tlist)+phase))])
        for i in range(1,order+1):
            c=0.5*hbar*omega0*(eta**i)*(1j)**i/math.factorial(i)
            c2=0.5*hbar*omega0*(eta**i)*(-1j)**i/math.factorial(i)
            H.append([c*tensor(s_plus,identity(2),destroy(size)**i),np.exp(1j*((detuning+i*nu)*tlist+phase))])
            H.append([c*tensor(identity(2),s_plus,destroy(size)**i),np.exp(1j*((detuning+i*nu)*tlist+phase))])
            H.append([c2*tensor(s_minus,identity(2),create(size)**i),np.exp(-1j*((detuning+i*nu)*tlist+phase))])
            H.append([c2*tensor(identity(2),s_minus,create(size)**i),np.exp(-1j*((detuning+i*nu)*tlist+phase))])
            H.append([c2*tensor(s_minus,identity(2),destroy(size)**i),np.exp(-1j*((detuning-i*nu)*tlist+phase))])
            H.append([c2*tensor(identity(2),s_minus,destroy(size)**i),np.exp(-1j*((detuning-i*nu)*tlist+phase))])
            H.append([c*tensor(s_plus,identity(2),create(size)**i),np.exp(1j*((detuning-i*nu)*tlist+phase))])
            H.append([c*tensor(identity(2),s_plus,create(size)**i),np.exp(1j*((detuning-i*nu)*tlist+phase))])
            for j in range(2,i):
                a=(destroy(size)**j)*create(size)**(i-j)
                H.append([c*tensor(s_plus,identity(2),a),np.exp(1j*detuning*tlist)])
                H.append([c*tensor(identity(2),s_plus,a),np.exp(1j*detuning*tlist)])
                H.append([c*tensor(s_minus,identity(2),a),np.exp(-1j*detuning*tlist)])
                H.append([c*tensor(identity(2),s_minus,a),np.exp(-1j*detuning*tlist)])
    if form ==1:# String format
        c=0.5*hbar*omega0
        H=[]
        H.append([c*tensor(s_plus,identity(size)),'exp(-1j*(detuning*t+phase))'])
        H.append([c*tensor(s_minus,identity(size)),'exp(1j*(detuning*t+phase))'])
        for i in range(1,order+1):
            c=0.5*hbar*omega0*(eta**i)/math.factorial(i)
            H.append([c*tensor(s_plus,destroy(size)**i),'exp(-1j*((detuning+%d*nu)*t+phase))'%i])
            H.append([c*tensor(s_minus,create(size)**i),'exp(1j*((detuning+%d*nu)*t+phase))'%i])
            H.append([c*tensor(s_minus,destroy(size)**i),'exp(1j*((detuning-%d*nu)*t+phase))'%i])
            H.append([c*tensor(s_plus,create(size)**i),'exp(-1j*((detuning-%d*nu)*t+phase))'%i])
            for j in range(2,i):
                a=(destroy(size)**j)*create(size)**(i-j)
                H.append([c*tensor(s_plus,a),'1'])
                H.append([c*tensor(s_minus,a),'1'])
    return H

def H_sideband(tlist,nu,omega0,eta,detuning,phase,size,atoms=[1,0],order=5,form=0):
    if form ==0:# Array format
        a=destroy(size)
        ad=create(size)
        expansion=[[],#0
                   [],#1
                   [[],[a*ad+ad*a]],#2
                   [[],[a*a*ad+a*ad*a+ad*a*a],
                    [a*ad*ad+ad*a*ad+ad*ad*a]],#3
                   [[],[a*a*a*ad+a*a*ad*a+a*ad*a*a+ad*a*a*a], 
                    [a*a*ad*ad+a*ad*a*ad+ad*a*a*ad+a*ad*ad*a+ad*a*ad*a+ad*ad*a*a],
                    [a*ad*ad*ad+ad*a*ad*ad+ad*ad*a*ad+ad*ad*ad*a]],#4
                   [[],[a*a*a*a*ad+a*a*a*ad*a+a*a*ad*a*a+a*ad*a*a*a+ad*a*a*a*a],
                    [a*a*a*ad*ad+a*a*ad*a*ad+a*ad*a*a*ad+ad*a*a*a*ad+a*a*ad*ad*a+a*ad*a*ad*a+ad*a*a*ad*a+a*ad*ad*a*a+ad*a*ad*a*a+ad*ad*a*a*a],
                    [a*a*ad*ad*ad+a*ad*a*ad*ad+ad*a*a*ad*ad+a*ad*ad*a*ad+ad*a*ad*a*ad+ad*ad*a*a*ad+a*ad*ad*ad*a+ad*a*ad*ad*a+ad*ad*a*ad*a+ad*ad*ad*a*a],
                    [a*ad*ad*ad*ad+ad*a*ad*ad*ad+ad*ad*a*ad*ad+ad*ad*ad*a*ad+ad*ad*ad*ad*a]],#5
                    ]# cross terms in taylor expansion
        sp=identity(2)# sigma_plus starting from identiy operation on 1st ion
        if atoms[1] == 0:
            sp=s_plus # sigma_plus starting from 1st ion
        a=identity(2)
        for i in range(1,atoms[0]):# generate series of identity operations on atoms other than desired one
            if i == atoms[1]:
                sp=tensor(sp,s_plus)
            else:
                sp=tensor(sp,identity(2))
            a=tensor(a,identity(2))

        a=tensor(a,destroy(size)) # annihilate operation, with n identity operation on atomic states
        c=0.5*hbar*omega0
        H=[]
        H.append([c*tensor(sp,identity(size)),np.exp(1j*((detuning*tlist)+phase))])
        H.append([c*tensor(sp.dag(),identity(size)),np.exp(-1j*((detuning*tlist)+phase))])
        for i in range(1,order+1):
            blue=tensor(sp,create(size)**i) # blue sideband
            red=tensor(sp,destroy(size)**i)
            
            c=0.5*hbar*omega0*((1j*eta)**i)/math.factorial(i)
            c2=0.5*hbar*omega0*((-1j*eta)**i)/math.factorial(i)
            H.append([red,c*np.exp(-1j*((detuning+i*nu)*tlist+phase))])
            H.append([red.dag(),c2*np.exp(1j*((detuning+i*nu)*tlist+phase))])
            H.append([blue.dag(),c2*np.exp(1j*((detuning-i*nu)*tlist+phase))])
            H.append([blue,c*np.exp(-1j*((detuning-i*nu)*tlist+phase))])
            for j in range(1,i):
                n=expansion[i][j][0]
                H.append([tensor(sp,n),c*np.exp(1j*(detuning+(j-i+j)*nu)*tlist)])
                H.append([tensor(sp.dag(),n),c2*np.exp(-1j*(detuning-(2*j-i)*nu)*tlist)])
    return H

def Multitone(tlist,nu,omega,eta,detuning,phase,size,order=5):
    H=[]
    for i in range(len(detuning)):
        H += H_sideband(tlist,nu,omega[i],eta,detuning[i],phase[i],size,[2,0],order)
        H += H_sideband(tlist,nu,omega[i],eta,-detuning[i],phase[i],size,[2,0],order)
        H += H_sideband(tlist,nu,omega[i],eta,detuning[i],phase[i],size,[2,1],order)
        H += H_sideband(tlist,nu,omega[i],eta,-detuning[i],phase[i],size,[2,1],order)
    return H

def MultitoneLD(tlist,nu,omega,eta,detuning,phase,size):
    jx=tensor(sigmax()/2,identity(2),identity(size))+tensor(identity(2),sigmax()/2,identity(size))
    jy=tensor(sigmay()/2,identity(2),identity(size))+tensor(identity(2),sigmay()/2,identity(size))
    x=tensor(identity(2),identity(2),(destroy(size)+create(size))/np.sqrt(2))
    p=tensor(identity(2),identity(2),1j*(-destroy(size)+create(size))/np.sqrt(2))
    
    H=[]
    H.append([2*omega0*jx,np.cos(detuning[0]*tlist+phase[0])])   
    for i in range(len(detuning)):
        H.append([2*omega[i]*jx,np.cos(detuning[i]*tlist+phase[i])])        
        H.append([-np.sqrt(2)*eta*omega[i]*jy*x,np.cos((nu-detuning[i])*tlist+phase[i])+np.cos((nu+detuning[i])*tlist+phase[i])])
        H.append([-np.sqrt(2)*eta*omega[i]*jy*p,np.sin((nu-detuning[i])*tlist+phase[i])+np.sin((nu+detuning[i])*tlist+phase[i])])
    return H
#%%
# MS Gate
order=1
atoms=2
hbar=1
size=10
n=2
nu=2*np.pi*200*10**3
omega0=nu*0.177
omega=np.array([1,-1,1,-1])*np.sqrt(3)*2/5*omega0
eta=0.1
detuning=nu*0.95
detunings=nu*(1-np.arange(0.05,0.25,0.05))
phase=np.zeros(4)
psi=basis(2,1)
for i in range(1,atoms):
    psi=tensor(psi,basis(2,1))
psi=tensor(psi,basis(size,2))


trange=np.pi/omega0*100
tlist=np.linspace(0,trange,int(100*trange*omega0/np.pi))#np.arange(0,trange,np.pi/omega0/100)
args={'detuning':detuning,'nu':nu,'phase':phase}
rho=[]
rhoj=[]
rhon=[]
H=Multitone(tlist,nu,omega,eta,detunings,phase,size,order)
#H=[]
#for i in range(atoms):
#    H=H+(H_sideband(tlist,nu,omega0,eta,detuning,phase,size,atoms=[atoms,i],order=order,form=0)
#        +H_sideband(tlist,nu,omega0,eta,-detuning,phase,size,atoms=[atoms,i],order=order,form=0))

sol=mesolve(H,psi,tlist,[],args=args,options=options)
for i in range(len(tlist)):
    rho.append(np.real(sol.states[i].ptrace([0,1])))
    rhoj.append(np.imag(sol.states[i].ptrace([0,1])))
    #rhon.append(sol.states[i].ptrace(2))
rho=np.array(rho)
rhoj=np.array(rhoj)



plt.figure()
for i in range(2):
    for j in range(2):
        plt.plot(tlist*nu,rho[:,i*2+j,i*2+j],label='|%s,%s>'%(['e','g'][i],['e','g'][j]))
plt.xlabel('Time (nu*t)')
plt.ylabel('Population')
plt.title('nu=200 kHz*2pi, omega0=0.177*nu, eta=0.1, detuing=0.95*nu+0.9*nu initial=|gg,2>')
plt.legend()

plt.figure(figsize=[8,8])
xev=np.linspace(-5,5,200)
plt.contourf(xev,xev,wigner(sol.states[-1].ptrace(2),xev,xev),100)
plt.xlim([-5,5])
plt.ylim([-5,5])
plt.xlabel('x')
plt.ylabel('p')
plt.title('final motional state')
#plt.figure()
#plt.plot(tlist*nu,states[:,0,0,2],label='|g,g,2>')
#plt.plot(tlist*nu,states[:,0,0,3],label='|g,g,3>')
#plt.plot(tlist*nu,states[:,1,0,3],label='|e,g,3>')
#plt.plot(tlist*nu,states[:,0,1,3],label='|g,e,3>')
#plt.plot(tlist*nu,states[:,1,1,2],label='|e,e,2>')


#%%
#fast gate
order=1
hbar=1
size=10
n=np.sqrt(2)
nu=2*np.pi*200*10**3
omega0=nu*0.177
eta=0.1
detuning=nu*0.95
psi=tensor(basis(2,1),basis(2,1),basis(size,2))#Qobj(thermal_dm(size,2)*(np.zeros([size,1])+1)))
phase=0
trange=np.pi/omega0*100#np.pi*(nu-detuning)/2/(eta*omega0)**2
tlist=np.linspace(0,trange,int(100*trange*omega0/np.pi))#np.arange(0,trange,np.pi/omega0/100)
args={'detuning':detuning,'nu':nu,'phase':phase}

jx=tensor(sigmax()/2,identity(2),identity(size))+tensor(identity(2),sigmax()/2,identity(size))
jy=tensor(sigmay()/2,identity(2),identity(size))+tensor(identity(2),sigmay()/2,identity(size))
x=tensor(identity(2),identity(2),(destroy(size)+create(size))/np.sqrt(2))
p=tensor(identity(2),identity(2),1j*(-destroy(size)+create(size))/np.sqrt(2))
H=[]
H.append([2*np.sqrt(2/3)*omega0*jx,np.cos(detuning*tlist)])
H.append([-np.sqrt(4/3)*eta*omega0*jy*x,np.cos((nu-detuning)*tlist)+np.cos((nu+detuning)*tlist)])
H.append([-np.sqrt(4/3)*eta*omega0*jy*p,np.sin((nu-detuning)*tlist)+np.sin((nu+detuning)*tlist)])
H.append([np.sqrt(4/3)*eta*omega0*jy*x,np.cos(2*(nu-detuning)*tlist)+np.cos(2*(nu+detuning)*tlist)])
H.append([np.sqrt(4/3)*eta*omega0*jy*p,np.sin(2*(nu-detuning)*tlist)+np.sin(2*(nu+detuning)*tlist)])
#H.append([-np.sqrt(2)*omega0*eta*jy*x,np.cos((nu-detuning)*tlist)])
#H.append([-np.sqrt(2)*omega0*eta*jy*p,np.sin((nu-detuning)*tlist)])

rho=[]
rhoj=[]
F=1-2*omega0**2/2/detuning**2*(1-np.cos(2*detuning*tlist)) # Fidelity

sol=mesolve(H,psi,tlist,[],[],args=args,options=options)
for i in range(len(tlist)):
    rho.append(np.real(sol.states[i].ptrace([0,1])))
    rhoj.append(np.imag(sol.states[i].ptrace([0,1])))

rho=np.array(rho)
rhoj=np.array(rhoj)

plt.figure()
for i in range(2):
    for j in range(2):
        plt.plot(tlist*nu,rho[:,i*2+j,i*2+j],label='|%s,%s>'%(['e','g'][i],['e','g'][j]))
plt.xlabel('Time (nu*t)')
plt.ylabel('Population')
plt.title('nu=200 kHz*2pi, omega0=0.177*nu, eta=0.1, detuing=0.95*nu+0.9*nu initial=|gg,2>')
plt.legend()

plt.figure(figsize=[8,8])
xev=np.linspace(-5,5,200)
plt.contourf(xev,xev,wigner(sol.states[-1].ptrace(2),xev,xev),100)
plt.xlim([-5,5])
plt.ylim([-5,5])
plt.xlabel('x')
plt.ylabel('p')
plt.title('final motional state')
#%%
w_gg=[]
w_eg=[]
w_ge=[]
w_ee=[]
frame_len=int(np.pi*(nu-detuning)/2/(eta*omega0)**2/tlist[1]*2)
for i in np.arange(0,frame_len,int(frame_len/10)):
    print(i)
    w_gg.append(wigner(tensor(basis(2,1),basis(2,1),qeye(size)).dag()*sol.states[i],xev,xev))
    w_eg.append(wigner(tensor(basis(2,0),basis(2,1),qeye(size)).dag()*sol.states[i],xev,xev))
    w_ge.append(wigner(tensor(basis(2,1),basis(2,0),qeye(size)).dag()*sol.states[i],xev,xev))
    w_ee.append(wigner(tensor(basis(2,0),basis(2,0),qeye(size)).dag()*sol.states[i],xev,xev))
#%%
rr=tensor(basis(2,1)+1j*basis(2,0),basis(2,1)+1j*basis(2,0),qeye(size))/2
ll=tensor(basis(2,1)-1j*basis(2,0),basis(2,1)-1j*basis(2,0),qeye(size))/2
rl=tensor(basis(2,1)+1j*basis(2,0),basis(2,1)-1j*basis(2,0),qeye(size))/2
lr=tensor(basis(2,1)-1j*basis(2,0),basis(2,1)+1j*basis(2,0),qeye(size))/2
op=[rr,rl,lr,ll]
for i in range(2):
    for j in range(2):
        frame_len=int(np.pi*(nu-detuning)/2/(eta*omega0)**2/tlist[1])
        # Animation of Phase Space
        fig=plt.figure(figsize=[8,8])
        ax = fig.add_subplot(1,1,1)
        
        # Setting limits for x and y axis
        ax.set_xlim(-5, 5)
        ax.set_ylim(-5, 5)
        
        def animation_function(frame):
            #where=np.where(w_gg[frame]==np.max(w_gg[frame]))
            #plt.scatter(xev[where[0]],xev[where[1]])
            plt.contourf(xev,xev,wigner(op[2*i+j].dag()*sol.states[frame],xev,xev),100)
            #plt.contourf(xev,xev,wigner(sol.states[frame].ptrace(2),xev,xev),100)
          
        animation = FuncAnimation(fig,
                                  func = animation_function,
                                  frames = np.arange(0, frame_len,int(frame_len/20)), 
                                  interval = 100 # Delay between frames in milliseconds.
                                  )
        plt.title('Motional State for |%s%s>'%(['+','-'][i],['+','-'][j]))
        plt.xlabel('x')
        plt.ylabel('p')
        animation.save('4tones_%s%s.gif'%(['right','left'][i],['right','left'][j]))
        plt.show()


#%%
# varying eta 
order=5
atoms=2
hbar=1
size=10
n=2
nu=2*np.pi*200*10**3
omega0=nu*0.1
eta=0.1
detuning=nu*0.95

psi=basis(2,1)
for i in range(1,atoms):
    psi=tensor(psi,basis(2,1))
psi=tensor(psi,basis(size,n))

phase=0
trange=np.pi/omega0*100
tlist=np.linspace(0,trange,int(100*trange*omega0/np.pi))#np.arange(0,trange,np.pi/omega0/100)
args={'detuning':detuning,'nu':nu,'phase':phase}
for eta in np.linspace(0.1,0.3,11):
    print(eta)
    rho=[]
    rhoj=[]
    rhon=[]
    
    H=[]
    #H=H_sideband2(tlist,nu,omega0,eta,detuning,phase,size,order,form=0)+H_sideband2(tlist,nu,omega0,eta,-detuning,phase,size,order,form=0)
    for i in range(atoms):
        H=H+(H_sideband(tlist,nu,omega0,eta,detuning,phase,size,atoms=[atoms,i],order=order,form=0)
            +H_sideband(tlist,nu,omega0,eta,-detuning,phase,size,atoms=[atoms,i],order=order,form=0))
    
    sol=mesolve(H,psi,tlist,[],args=args,options=options)
    for i in range(len(tlist)):
        rho.append(np.real(sol.states[i].ptrace([0,1])))
        rhoj.append(np.imag(sol.states[i].ptrace([0,1])))
        #rhon.append(sol.states[i].ptrace(2))
    rho=np.array(rho)
    rhoj=np.array(rhoj)
    
    
    
    plt.figure()
    for i in range(2):
        for j in range(2):
            plt.plot(tlist*nu,rho[:,i*2+j,i*2+j],label='|%s,%s>'%(['e','g'][i],['e','g'][j]))
    plt.xlabel('Time (nu*t)')
    plt.ylabel('Population')
    plt.title('nu=200 kHz*2pi, omega0=0.177*nu, eta=%s, detuing=0.95*nu initial=|gg,2>'%(eta))
    plt.legend()
    plt.savefig('MS_200k_0100_0%d_095_gg2.jpg'%(eta*100))
    # name=('trap_frequency(2*pi), omega0(nu), eta, detuning(nu), size, initial state')


