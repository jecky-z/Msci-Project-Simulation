# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 21:14:30 2021

@author: jecky
"""

#%%
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import math
from qutip import *
import timeit
options=qutip.solver.Options(atol=1e-07,rtol=1e-07,store_states=True,num_cpus=8)
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
hbar=1
s_plus=(sigmax()+1j*sigmay())*0.5
s_minus=(sigmax()-1j*sigmay())*0.5

def H_sideband(nu,omega0,eta,detuning,size,order=5,form=0):
    if form ==0:# Array format
        c=0.5*hbar*omega0
        H=[]
        H.append([c*tensor(s_plus,qeye(size)),np.exp(-1j*(detuning)*tlist)])
        H.append([c*tensor(s_minus,qeye(size)),np.exp(1j*(detuning)*tlist)])
        for i in range(1,order+1):
            c=0.5*hbar*omega0*(eta**i)/math.factorial(i)
            H.append([c*tensor(s_plus,destroy(size)**i),np.exp(-1j*(detuning+i*nu)*tlist)])
            H.append([c*tensor(s_minus,create(size)**i),np.exp(1j*(detuning+i*nu)*tlist)])
            H.append([c*tensor(s_minus,destroy(size)**i),np.exp(1j*(detuning-i*nu)*tlist)])
            H.append([c*tensor(s_plus,create(size)**i),np.exp(-1j*(detuning-i*nu)*tlist)])
            for j in range(2,i):
                a=(destroy(size)**j)*create(size)**(i-j)
                H.append([c*tensor(s_plus,a),1+0*tlist])
                H.append([c*tensor(s_minus,a),1+0*tlist])
    if form ==1:# String format
        c=0.5*hbar*omega0
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
    return H
#%%
order=5
hbar=1
size=10
n=2
nu=2*np.pi*1.1*10**6
omega0=nu/5#2*np.pi*6.94*10**3
eta=0.09
detuning=nu*1
psi=tensor(basis(2,1),basis(size,n))
trange=np.pi/omega0*1
tlist=np.arange(0,trange,np.pi/omega0/100)

shift=np.zeros([50,50])
for o in range(1,50):
    for delt in range(20,21):
        print([o,delt])
        n_scan=100
        scan=(np.arange(0,1,1/n_scan)-0.5)
        c=[]
        r=(omega0*(o*0.02))**2/(2*delt*0.1)/nu**2
        for i in range(n_scan):
            #print(i)
            detuning=scan[i]*nu*4*r
            H=H_sideband(nu,omega0,eta,detuning,size,order,form=0)+H_sideband(nu,omega0*(o*0.02),eta,nu*(delt*0.1),size,order,form=0)
            sol=mesolve(H,psi,tlist,[],[tensor(sigmaz(),qeye(size)),tensor(qeye(2),num(size))],args={'detuning':detuning,'nu':nu},options=options)
            c.append(sol)
        
        expect=np.zeros([n_scan,2,len(tlist)])
        states=np.zeros([n_scan,len(tlist),2*size,1])
        for i in range(n_scan):
            expect[i]=c[i].expect
            for j in range(len(tlist)):
                states[i,j]=c[i].states[j].full()*(np.conj(c[i].states[j].full())) # extract individual states
        shift[o,delt]=scan[np.where(expect[:,0,-1]==max(expect[:,0,-1]))][0]*4*r
#np.savetxt('0to5_0to5.txt',shift)
plt.figure()
plt.plot(scan*0.4,expect[:,0,-1])
plt.xlabel('detuing/nu')
plt.ylabel('Expectation sigma_z')
plt.title('nu=%d *2pi, omega0=nu/%d, eta=%s, initial=|g,%d>'%(round(nu/2/np.pi),round(nu/omega0),eta,n))
#plt.figure()
#plt.plot(scan*1,expect[:,1,-1])
#plt.xlabel('detuing/nu')
#plt.title('Expectation n')
#%%
order=5
hbar=1
size=10
n=2
nu=2*np.pi*1.1*10**6
omega0=nu/5#2*np.pi*6.94*10**3
omega2=omega0*2
detuning2=nu*4
eta=0.09
detuning=nu*1
psi=basis(2,1)
trange=np.pi/omega0*1
tlist=np.arange(0,trange,np.pi/omega0/100)


n_scan=1000
scan=(np.arange(0,1,1/n_scan)-0.5)
result=[]
for i in range(n_scan):
    print(i)
    detuning=scan[i]*nu*1
    c=0.5*hbar*omega0
    d=0.5*hbar*omega2
    H=[]
    H.append([c*tensor(s_plus,qeye(size)),np.exp(-1j*(detuning)*tlist)])
    H.append([d*tensor(s_plus,qeye(size)),np.exp(-1j*(detuning2)*tlist)])
    H.append([c*tensor(s_minus,qeye(size)),np.exp(1j*(detuning)*tlist)])
    H.append([d*tensor(s_minus,qeye(size)),np.exp(1j*(detuning2)*tlist)])

    sol=mesolve(H,psi,tlist,[],[tensor(sigmaz(),qeye(size)),tensor(qeye(2),num(size))],args={'detuning':detuning,'nu':nu},options=options)
    result.append(sol)

expect=np.zeros([n_scan,2,len(tlist)])
states=np.zeros([n_scan,len(tlist),2*size,1])
for i in range(n_scan):
    expect[i]=result[i].expect
    for j in range(len(tlist)):
        states[i,j]=result[i].states[j].full()*(np.conj(result[i].states[j].full())) # extract individual states

plt.figure()
plt.plot(scan*1,expect[:,0,-1])
plt.xlabel('detuing/nu')
plt.ylabel('Expectation sigma_z')
plt.title('nu=%d *2pi, omega0=nu/%d, eta=%s, initial=|g,%d>'%(round(nu/2/np.pi),round(nu/omega0),eta,n))


















