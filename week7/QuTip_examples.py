# -*- coding: utf-8 -*-
"""
Created on Sat Nov  6 19:45:33 2021

@author: jecky
"""

#%%
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import math
from qutip import *
import timeit
options=qutip.solver.Options(atol=1e-07,rtol=1e-07,store_states=True)
params = {
   'axes.labelsize': 25,
   'font.size': 25,
   'legend.fontsize': 18,
   'xtick.labelsize': 25,
   'ytick.labelsize': 25,
   'figure.figsize': [20,20]
   } 
plt.rcParams.update(params)
name=('trap_frequency(2*pi), omega0(nu), eta, detuning(nu), size, initial state')
s_plus=(sigmax()+1j*sigmay())*0.5
s_minus=(sigmax()-1j*sigmay())*0.5

def H_sideband(nu,omega0,eta,detuning,phase=0,size=10,order=5,form=0):
    if form ==0:# Array format
        c=0.5*hbar*omega0
        H=[]
        H.append([c*tensor(s_plus,qeye(size)),np.exp(-1j*((detuning*tlist)+phase))])
        H.append([c*tensor(s_minus,qeye(size)),np.exp(1j*((detuning*tlist)+phase))])
        for i in range(1,order+1):
            c=0.5*hbar*omega0*(eta**i)/math.factorial(i)
            H.append([c*tensor(s_plus,destroy(size)**i),np.exp(-1j*((detuning+i*nu)*tlist+phase))])
            H.append([c*tensor(s_minus,create(size)**i),np.exp(1j*((detuning+i*nu)*tlist+phase))])
            H.append([c*tensor(s_minus,destroy(size)**i),np.exp(1j*((detuning-i*nu)*tlist+phase))])
            H.append([c*tensor(s_plus,create(size)**i),np.exp(-1j*((detuning-i*nu)*tlist+phase))])
            for j in range(2,i):
                a=(destroy(size)**j)*create(size)**(i-j)
                H.append([c*tensor(s_plus,a),1+0*tlist])
                H.append([c*tensor(s_minus,a),1+0*tlist])
    if form ==1:# String format
        c=0.5*hbar*omega0
        H=[]
        H.append([c*tensor(s_plus,qeye(size)),'exp(-1j*(detuning*t+phase))'])
        H.append([c*tensor(s_minus,qeye(size)),'exp(1j*(detuning*t+phase))'])
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
#%%
# Sideband Scan
hbar=1
size=10
n=2
nu=2*np.pi*1.1*10**6 # 1.1 MHz
omega0=2*np.pi*6.94*10**3/0.09
eta=0.1
detuning=nu*1
phase=0
psi=tensor(basis(2,1),basis(size,n))# used (0,1) as ground state, eigenvalue sigma_z of -1
trange=np.pi/omega0*1
tlist=np.arange(0,trange,trange/100)
args={'detuning':detuning,'nu':nu,'phase':phase}
n_scan=100
scan=(np.arange(0,1,1/n_scan)+1)
c=[]
for i in range(n_scan):
    print(i)
    detuning=scan[i]*nu*1
    H=H_sideband(nu,omega0,eta,detuning,phase,size,order=5,form=0)#+H_sideband(nu,omega0,eta,nu*-25,size,order=5,form=0)
    sol=mesolve(H,psi,tlist,[],[tensor(sigmaz(),qeye(size)),tensor(qeye(2),num(size))],args=args,options=options)
    c.append(sol)


expect=np.zeros([n_scan,2,len(tlist)])
states=np.zeros([n_scan,len(tlist),2*size,1])
for i in range(n_scan):
    expect[i]=c[i].expect
    for j in range(n_scan):
        states[i,j]=c[i].states[j].full()*(np.conj(c[i].states[j].full())) # extract individual states
plt.figure()
plt.plot(scan*1,expect[:,0,-1])
plt.xlabel('detuing/nu')
plt.title('Expectation sigma_z')
plt.figure()
plt.plot(scan*1,expect[:,1,-1])
plt.xlabel('detuing/nu')
plt.title('Expectation n')
#%%
# plot individual states
plt.figure()
for i in range(2):
    for j in range(size):
        plt.plot(tlist*omega0/np.pi,states[50,:,i*10+j],label='|%s,%s>'%(i,j))
plt.legend()
plt.xlabel('Time (pi/omega0)')
plt.title('population vs t')
print('run time=',start-end,'s')
#%%
# Sideband transition using array format
order=5
hbar=1
size=10
n=2
nu=2*np.pi*1.1*10**6
omega0=2*np.pi*6.94*10**3/0.09
eta=0.09
detuning=nu*1
psi=tensor(basis(2,1),basis(size,n))#+tensor(basis(2,1),basis(size,1))+tensor(basis(2,1),basis(size,0)) # used (0,1) as ground state, eigenvalue sigma_z of -1
phase=0
trange=np.pi/omega0*100
tlist=np.arange(0,trange,np.pi/omega0/100)
s_plus=(sigmax()+1j*sigmay())*0.5
s_minus=(sigmax()-1j*sigmay())*0.5
args={'detuning':detuning,'nu':nu,'phase':phase}
H=H_sideband(nu,omega0,eta,detuning,phase,size,order=5,form=0)#+H_sideband(nu,omega0*5,eta,nu*0,size,order=0,form=0)

sol=mesolve(H,psi,tlist,[],[tensor(sigmaz(),qeye(size)),tensor(qeye(2),num(size))],args=args,options=options)

plt.figure()
plt.plot(tlist*omega0/np.pi,sol.expect[0])
plt.xlabel('Time (pi/omega)')
plt.ylabel('Expectation sigma_z')
plt.title('nu=%d *2pi, omega0=nu/%d, eta=%s, initial=|g,%d>'%(round(nu/2/np.pi),round(nu/omega0),eta,n))
#%%
# Using string format
order=5
hbar=1
size=10
n=2
nu=2*np.pi*10**3
omega0=nu/5
eta=0.1
detuning=nu*2
psi=tensor(basis(2,1),basis(size,n)) # used (0,1) as ground state, eigenvalue sigma_z of -1
s_plus=(sigmax()+1j*sigmay())*0.5
s_minus=(sigmax()-1j*sigmay())*0.5
trange=np.pi/omega0*100
tlist=np.arange(0,trange,np.pi/omega0/100)
s_plus=(sigmax()+1j*sigmay())*0.5
s_minus=(sigmax()-1j*sigmay())*0.5

c=0.5*hbar*omega0
args={'detuning':detuning,'nu':nu}

H=[]
H.append([c*tensor(s_plus,qeye(size)),'exp(-1j*detuning*t)'])
H.append([c*tensor(s_minus,qeye(size)),'exp(1j*detuning*t)'])
for i in range(1,order+1):
    c=0.5*hbar*omega0*(eta**i)/math.factorial(i)
    H.append([c*tensor(s_plus,destroy(size)**i),'exp(-1j*(detuning+%d*nu)*t)'%i])
    H.append([c*tensor(s_minus,create(size)**i),'exp(1j*(detuning+%d*nu)*t)'%i])
    H.append([c*tensor(s_minus,destroy(size)**i),'exp(1j*(detuning-%d*nu)*t)'%i])
    H.append([c*tensor(s_plus,create(size)**i),'exp(-1j*(detuning-%d*nu)*t)'%i])
    for j in range(2,i):
        a=(destroy(size)**j)*create(size)**(i-j)
        H.append([c*tensor(s_plus,a),'1'])
        H.append([c*tensor(s_minus,a),'1'])
        
start=timeit.default_timer()
sol=mesolve(H,psi,tlist,[],[tensor(sigmaz(),qeye(size)),tensor(qeye(2),num(size))],args=args,options=options)
end=timeit.default_timer()
plt.figure()
plt.plot(tlist*omega0/np.pi,sol.expect[0])
plt.xlabel('Time (pi/omega)')
plt.ylabel('Expectation sigma_z')
plt.title('nu=%d *2pi, omega0=nu/%d, eta=%s, initial=|g,%d>'%(round(nu/2/np.pi),round(nu/omega0),eta,n))
#%%
# A sequency of lasers
order=5
hbar=1
size=10
n=0
nu=2*np.pi*1.1*10**6
omega0=2*np.pi*6.94*10**3/0.09
eta=0.09
omega1=omega0*eta
detuning=[0,-1,0,-1,-1,0,-1,0,-1]
psi=tensor(basis(2,1),basis(size,0))#+tensor(basis(2,1),basis(size,1)))/np.sqrt(2)
phase=[0,-0.5,1,0.5,0,0.71,-0.29,0.1,-0.51]
trange=[np.pi/omega0*0.5,np.pi/omega1*0.7,np.pi/omega0*0.73,np.pi/omega1*0.71,
        np.pi/omega1*0.71,np.pi/omega0*0.5,np.pi/omega1*1.42,np.pi/omega0*1.59,np.pi/omega1*0.72]
expect=[]
t=[[0]]
for i in range(9):
    print(i)
    tlist=np.arange(0,trange[i],np.pi/omega0/100)
    s_plus=(sigmax()+1j*sigmay())*0.5
    s_minus=(sigmax()-1j*sigmay())*0.5
    args={'detuning':detuning[i]*nu,'nu':nu,'phase':phase}
    H=H_sideband(nu,omega0,eta,detuning[i]*nu,phase[i],size,order,form=0)#+H_sideband(nu,omega0*5,eta,nu*0,size,order=0,form=0)
    
    sol=mesolve(H,psi,tlist,[],[tensor(sigmaz(),qeye(size))],args=args,options=options)
    expect.append(sol.expect)
    t.append(tlist+t[-1][-1])
    psi=sol.states[-1]

expect=np.array(expect).reshape([-1])
t=np.array(t)
plt.figure()
for i in range(9):
    plt.plot(t[i+1],0.5*(expect[i]+1))
plt.xlabel('Time (s)')
plt.ylabel('Expectation sigma_z')
plt.title('sequence of pulses,nu=1.1MHz*2pi, omega0=77kHz*2pi, eta=0.09, initial=|g,0>')  
    
    