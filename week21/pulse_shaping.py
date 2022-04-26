# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 16:02:31 2022

@author: jecky
"""

#%%
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import math
from qutip import *
import itertools

options=qutip.solver.Options(atol=1e-12,rtol=1e-12,store_states=True)
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


def H_sideband(tlist,nu,omega0,eta,detuning,phase,size,atoms=[1,0],order=5,envelope=1):
    a=destroy(size)
    ad=create(size)
    expansion=[]
    for i in range(order+1): # the order of expansion
        sub=[]
        sub.append([])
        for j in range(1,i): # the numer of ad in the cross terms
            op=[ad]*j+[a]*(i-j)
            cross_terms_list=list(itertools.permutations(op))
            cross_term=0
            for k in range(len(cross_terms_list)):
                dummy=1
                for l in range(i):
                    dummy *= cross_terms_list[k][l]
                cross_term += dummy
            sub.append([cross_term/math.factorial(j)/math.factorial(i-j)])
        expansion.append(sub)

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
    c=0.5*hbar*omega0*envelope
    H=[]
    H.append([tensor(sp,identity(size)),c*np.exp(-1j*((detuning*tlist)+phase))])
    H.append([tensor(sp.dag(),identity(size)),c*np.exp(1j*((detuning*tlist)+phase))])
    for i in range(1,order+1):
        blue=tensor(sp,create(size)**i) # blue sideband
        red=tensor(sp,destroy(size)**i)
        
        c=0.5*hbar*omega0*((1j*eta)**i)/math.factorial(i)*envelope
        c2=0.5*hbar*omega0*((-1j*eta)**i)/math.factorial(i)*envelope
        H.append([red,c*np.exp(-1j*((detuning+i*nu)*tlist+phase))])
        H.append([red.dag(),c2*np.exp(1j*((detuning+i*nu)*tlist+phase))])
        H.append([blue.dag(),c2*np.exp(1j*((detuning-i*nu)*tlist+phase))])
        H.append([blue,c*np.exp(-1j*((detuning-i*nu)*tlist+phase))])
        for j in range(1,i):
            n=expansion[i][j][0]
            H.append([tensor(sp,n),c*np.exp(-1j*(detuning+(j-i+j)*nu)*tlist)])
            H.append([tensor(sp.dag(),n),c2*np.exp(1j*(detuning-(2*j-i)*nu)*tlist)])
    return H
#%%
# Multi Sideband 
order=5
atoms=2
hbar=1
size=10
n=0
nu=2*np.pi*461*10**3
eta=0.1
omega0=nu*0.079
omega=omega0

detuning=0.98*nu

phase=0
psi=tensor(basis(2,1),basis(2,1),Qobj(thermal_dm(size,n)*(np.zeros([size,1])+1)))
#tensor(ms.Sy.eigenstates()[1][0],basis(size,0))
#basis(2,1).proj(),basis(2,1).proj()

tw=2*np.pi/nu*5
trange=np.pi*(2)/delta
tstep=np.pi/delta/4000
tlist=np.arange(0,trange+tw,tstep)
turnon=np.where(tlist<=tw/2+tstep)[0][-1]
tp=np.where(tlist<=trange+tw/2)[0][-1]
gate_time=np.where(tlist <= trange+tw/2)[0][-1]
args={'detuning':detuning,'nu':nu,'phase':phase}
envelope=np.concatenate([np.sin(np.pi*tlist[:turnon]/tw)**2,
                   1+0*tlist[turnon:tp],
                   np.sin(np.pi*(tlist[tp:]-trange)/tw)**2])
#envelope=np.ones(len(tlist))
H=H_sideband(tlist,nu,omega0,eta,detuning,phase,size,atoms=[2,0],order=order,envelope=envelope)





target=tensor(basis(2,0),basis(2,1))
target=tensor(target,qeye(size))
t_rho=target*target.dag()

sol=mesolve(H,psi,tlist,[],[t_rho],args=args,options=options)
#(rr*rr.dag())*x,(rr*rr.dag())*p

rho=[]
rhoj=[]
for i in range(len(tlist)):
    rho.append(np.real(sol.states[i].ptrace([0,1])))
    rhoj.append(np.imag(sol.states[i].ptrace([0,1])))
    #rhon.append(sol.states[i].ptrace(2))
rho=np.array(rho)
rhoj=np.array(rhoj)



plt.figure()
for i in range(2):
    for j in range(2):
        plt.plot((tlist-0.5*tw)/trange,rho[:,i*2+j,i*2+j],label='|%s,%s>'%(['e','g'][i],['e','g'][j]))
plt.plot((tlist-0.5*tw)/trange,envelope,label='envelope')
plt.xlabel('Gate Time')
plt.ylabel('Population')
plt.title('nu=%d kHz*2pi, omega0=%s*nu, eta=%s, n=%s'%(nu/2000/np.pi,omega0/nu,eta,n))
plt.legend()
#%%
# MS Gate
# Multi Sideband 
order=5
atoms=2
hbar=1
size=100
n=0
eta=0.1
nu=2*np.pi*461*10**3
omega0=nu*0.079*0.5*eta*np.exp(0.5*eta**2)
omega=np.array([1,1])*omega0
phase=0

target=(tensor(basis(2,1),basis(2,1)) - 1j*tensor(basis(2,0),basis(2,0)))/np.sqrt(2)
target=tensor(target,qeye(size))
t_rho=target*target.dag()
F=[]
num=0
options=qutip.solver.Options(atol=1e-12,rtol=1e-12,nsteps=1e6,store_states=True)
loop=np.linspace(0.01,20,100)


for e in loop:
    print(e)
    delta=2*eta*omega0
    tw=2*np.pi/nu*e
    trange=np.pi*(2)/delta
    tstep=np.pi/delta/4000
    tlist=np.arange(0,trange+tw,tstep)
    turnon=np.where(tlist<=tw/2+tstep)[0][-1]
    tp=np.where(tlist<=trange+tw/2)[0][-1]
    gate_time=np.where(tlist <= trange+tw/2)[0][-1]
    args={'detuning':detuning,'nu':nu,'phase':phase}
    envelope=np.concatenate([np.sin(np.pi*tlist[:turnon]/tw)**2,
                       1+0*tlist[turnon:tp],
                       np.sin(np.pi*(tlist[tp:]-trange)/tw)**2])
    
    H=H_sideband(tlist,nu,omega0,eta,0.5*nu,phase,size,atoms=[2,0],order=order,envelope=envelope)
    target=tensor(basis(2,0),basis(2,1))
    target=tensor(target,qeye(size))
    t_rho=target*target.dag()
    
    sol=mesolve(H,psi,tlist,[],[t_rho],args=args,options=options)
    F.append(max(sol.expect[0]))


F=np.reshape(np.array(F),[-1])
plt.figure()
plt.plot(loop,F)

plt.xlabel('turn on time')
plt.ylabel('Infidelity')
plt.title('nu=461 kHz*2pi, omega0=%s*nu, eta=%s,n=%s'%(omega0/nu,eta,n))
plt.legend()
    








