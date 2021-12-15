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

        a=tensor(a,destroy(size)) # annihilate operation, with n identity operation on aatomic states
        c=0.5*hbar*omega0
        H=[]
        H.append([c*tensor(sp,identity(size)),np.exp(1j*((detuning*tlist)+phase))])
        H.append([c*tensor(sp.dag(),identity(size)),np.exp(-1j*((detuning*tlist)+phase))])
        for i in range(1,order+1):
            blue=tensor(sp,create(size)**i) # blue sideband
            red=tensor(sp,destroy(size)**i)

            c=0.5*hbar*omega0*((1j*eta)**i)/math.factorial(i)
            c2=0.5*hbar*omega0*((-1j*eta)**i)/math.factorial(i)
            H.append([c*red,np.exp(-1j*((detuning+i*nu)*tlist+phase))])
            H.append([c2*red.dag(),np.exp(1j*((detuning+i*nu)*tlist+phase))])
            H.append([c2*blue.dag(),np.exp(1j*((detuning-i*nu)*tlist+phase))])
            H.append([c*blue,np.exp(-1j*((detuning-i*nu)*tlist+phase))])
            for j in range(2,i):
                n=(destroy(size)**j)*(create(size)**(i-j))
                H.append([c*tensor(sp,n),np.exp(1j*detuning*tlist)])
                H.append([c*tensor(sp.dag(),n),np.exp(-1j*detuning*tlist)])
    return H
#%%
# MS Gate
order=1
atoms=2
hbar=1
size=200
n=2
nu=2*np.pi*200*10**3
omega0=nu*0.177
eta=0.1
detuning=nu*0.95

psi=basis(2,1)
for i in range(1,atoms):
    psi=tensor(psi,basis(2,1))
psi=tensor(psi,coherent(size,np.sqrt(0)))

phase=0
trange=np.pi/omega0*200
tlist=np.linspace(0,trange,int(100*trange*omega0/np.pi))#np.arange(0,trange,np.pi/omega0/100)
args={'detuning':detuning,'nu':nu,'phase':phase}
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
plt.title('nu=200 kHz*2pi, omega0=0.177*nu, eta=0.1, detuing=0.95*nu initial=|gg,alpha=sqrt(2)>')
plt.legend()

plt.figure()
xev=np.linspace(-5,5,200)
plt.contourf(xev,xev,wigner(sol.states[-1].ptrace(2),xev,xev),100)
plt.xlabel('Time (nu*t)')
plt.ylabel('Population')
plt.title('nu=200 kHz*2pi, omega0=0.1*nu, eta=0.1, initial=|eg>')
#plt.figure()
#plt.plot(tlist*nu,states[:,0,0,2],label='|g,g,2>')
#plt.plot(tlist*nu,states[:,0,0,3],label='|g,g,3>')
#plt.plot(tlist*nu,states[:,1,0,3],label='|e,g,3>')
#plt.plot(tlist*nu,states[:,0,1,3],label='|g,e,3>')
#plt.plot(tlist*nu,states[:,1,1,2],label='|e,e,2>')
#%%
# Retr
order=5
hbar=1
size=200
n=np.sqrt(2)
nu=2*np.pi*200*10**3
omega0=nu*0.05
eta=0.1
detuning=nu*-0.9
omega1=(omega0*eta)**2/(nu-abs(detuning))
psi=tensor((basis(2,1)+basis(2,0))/np.sqrt(2),basis(2,1),coherent(size,n))

phase=0
trange=np.pi/omega1
tlist=np.linspace(0,trange,int(100*trange*omega0/np.pi))#np.arange(0,trange,np.pi/omega0/100)
args={'detuning':detuning,'nu':nu,'phase':phase}
rho=[]
rhoj=[]

H=H_sideband2(tlist,nu,omega0,eta,detuning,phase,size,order,form=0)#+H_sideband2(nu,omega0,eta,-detuning,phase,size,order,form=0)

sol=mesolve(H,psi,tlist,[],[t],args=args,options=options)
for i in range(len(tlist)):
    rho.append(np.real(sol.states[i].ptrace([0,1])))
    rhoj.append(np.imag(sol.states[i].ptrace([0,1])))


H=H_sideband2(tlist+tlist[-1],nu,omega0,eta,-detuning,phase,size,order,form=0)
sol=mesolve(H,sol.states[-1],tlist+tlist[-1],[],[tensor(sigmaz(),identity(2),identity(size)),tensor(identity(2),sigmaz(),identity(size)),tensor(identity(2),identity(2),num(size))],args=args,options=options)

for i in range(len(tlist)):
    rho.append(np.real(sol.states[i].ptrace([0,1])))
    rhoj.append(np.imag(sol.states[i].ptrace([0,1])))
rho=np.array(rho)
rhoj=np.array(rhoj)

tlist=np.reshape(np.array([tlist,tlist+tlist[-1]]),[-1])
plt.figure()
for i in range(2):
    for j in range(2):
        plt.plot(tlist*nu,rho[:,i*2+j,i*2+j],label='|%s,%s>'%(['e','g'][i],['e','g'][j]))
plt.plot(tlist*nu,rho[:,0,1],label='Real:rho(ee,eg)')
plt.plot(tlist*nu,rhoj[:,0,1],label='Imag:rho(ee,eg)')
plt.xlabel('Time (nu*t)')
plt.ylabel('Population')
plt.ylim([-0.6,0.6])
plt.title('nu=200 kHz*2pi, omega0=0.05*nu, eta=0.1, initial=(|gg>+|eg>)/sqrt(2)')
plt.legend()
#%%
# Monte Carlo Method
order=5
hbar=1
size=10
n=2
nu=2*np.pi*200*10**3
omega0=nu*0.1
eta=0.1
detuning=nu*0.9
Gamma=nu*2*10**-4 # relaxation factor
n_therm=2
psi=tensor(basis(2,1),basis(2,1),basis(size,n))

phase=0
trange=np.pi/omega0*500
tlist=np.linspace(0,trange,int(100*trange*omega0/np.pi))#np.arange(0,trange,np.pi/omega0/100)
args={'detuning':detuning,'nu':nu,'phase':phase}
rho=[]
rhoj=[]

H=H_sideband2(tlist,nu,omega0,eta,detuning,phase,size,order,form=0)+H_sideband2(tlist,nu,omega0,eta,-detuning,phase,size,order,form=0)
C1=np.sqrt(Gamma*(1+n_therm))*tensor(identity(2),identity(2),destroy(size))
C2=np.sqrt(Gamma*(n_therm))*tensor(identity(2),identity(2),create(size))

col_times=[]
col_which=[]
for k in range(5):
    sol=mcsolve(H,psi,tlist,[C1,C2],[tensor(sigmaz(),identity(2),identity(size)),tensor(identity(2),sigmaz(),identity(size)),tensor(identity(2),identity(2),num(size))],
                args=args,ntraj=2,options=options)
    for i in range(2):
        col_times.append(sol.col_times[i])
        col_which.append(sol.col_which[0])
        for j in range(len(tlist)):
            rho.append(np.real(sol.states[i,j].ptrace([0,1])))
            rhoj.append(np.imag(sol.states[i,j].ptrace([0,1])))
rho=np.reshape(np.array(rho),[10,-1,4,4])
rhoj=np.reshape(np.array(rhoj),[10,-1,4,4])
plt.figure()
plt.plot(tlist*nu,np.sum(rho,axis=0)[:,0,0]/10,label='Average')

plt.xlabel('Time (nu*t)')
plt.ylabel('Population')
plt.title('nu=200 kHz*2pi, omega0=0.1*nu, eta=0.1, initial=|gg0>>')
plt.legend()
#%%
# C-NOT Gate
order=10
hbar=1
size=200
n=np.sqrt(2)
nu=2*np.pi*200*10**3
omega0=nu*0.05
eta=0.1
detuning=nu*-0.9
omega1=(omega0*eta)**2/(nu-abs(detuning))
psi=tensor((basis(2,1)+basis(2,0))/np.sqrt(2),basis(2,1),coherent(size,n))

phase=0
trange=np.pi/omega1
tlist=np.linspace(0,trange,int(100*trange*omega0/np.pi))#np.arange(0,trange,np.pi/omega0/100)
args={'detuning':detuning,'nu':nu,'phase':phase}
rho=[]
rhoj=[]

H=H_sideband2(tlist,nu,omega0,eta,detuning,phase,size,order,form=0)#+H_sideband2(nu,omega0,eta,-detuning,phase,size,order,form=0)

sol=mesolve(H,psi,tlist,[],[tensor(sigmaz(),identity(2),identity(size)),tensor(identity(2),sigmaz(),identity(size)),tensor(identity(2),identity(2),num(size))],args=args,options=options)
for i in range(len(tlist)):
    rho.append(np.real(sol.states[i].ptrace([0,1])))
    rhoj.append(np.imag(sol.states[i].ptrace([0,1])))


H=H_sideband2(tlist+tlist[-1],nu,omega0,eta,-detuning,phase,size,order,form=0)
sol=mesolve(H,sol.states[-1],tlist+tlist[-1],[],[tensor(sigmaz(),identity(2),identity(size)),tensor(identity(2),sigmaz(),identity(size)),tensor(identity(2),identity(2),num(size))],args=args,options=options)

for i in range(len(tlist)):
    rho.append(np.real(sol.states[i].ptrace([0,1])))
    rhoj.append(np.imag(sol.states[i].ptrace([0,1])))
rho=np.array(rho)
rhoj=np.array(rhoj)

tlist=np.reshape(np.array([tlist,tlist+tlist[-1]]),[-1])
plt.figure()
for i in range(2):
    for j in range(2):
        plt.plot(tlist*nu,rho[:,i*2+j,i*2+j],label='|%s,%s>'%(['e','g'][i],['e','g'][j]))
plt.plot(tlist*nu,rho[:,0,1],label='Real:rho(ee,eg)')
plt.plot(tlist*nu,rhoj[:,0,1],label='Imag:rho(ee,eg)')
plt.xlabel('Time (nu*t)')
plt.ylabel('Population')
plt.ylim([-0.6,0.6])
plt.title('nu=200 kHz*2pi, omega0=0.05*nu, eta=0.1, initial=(|gg>+|eg>)/sqrt(2)')
plt.legend()

#%%
#fast gate
order=10
hbar=1
size=200
n=np.sqrt(2)
nu=2*np.pi*200*10**3
omega0=nu*0.177
eta=0.1
detuning=nu*0.95
psi=tensor(basis(2,1),basis(2,1),Qobj(thermal_dm(size,2)*(np.zeros([size,1])+1)))
phase=0
trange=np.pi/omega0*200#np.pi*(nu-detuning)/2/(eta*omega0)**2
tlist=np.linspace(0,trange,int(100*trange*omega0/np.pi))#np.arange(0,trange,np.pi/omega0/100)
args={'detuning':detuning,'nu':nu,'phase':phase}

jx=tensor(sigmax()/2,identity(2),identity(size))+tensor(identity(2),sigmax()/2,identity(size))
jy=tensor(sigmay()/2,identity(2),identity(size))+tensor(identity(2),sigmay()/2,identity(size))
x=tensor(identity(2),identity(2),(destroy(size)+create(size))/np.sqrt(2))
p=tensor(identity(2),identity(2),1j*(-destroy(size)+create(size))/np.sqrt(2))
H=[]
#H.append([2*omega0*jx,np.cos(detuning*tlist)])
#H.append([-np.sqrt(2)*eta*omega0*jy*x,np.cos((nu-detuning)*tlist)+np.cos((nu+detuning)*tlist)])
#H.append([-np.sqrt(2)*eta*omega0*jy*p,np.sin((nu-detuning)*tlist)+np.sin((nu+detuning)*tlist)])
H.append([-np.sqrt(2)*omega0*eta*jy*x,np.cos((nu-detuning)*tlist)])
H.append([-np.sqrt(2)*omega0*eta*jy*p,np.sin((nu-detuning)*tlist)])

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
plt.title('nu=200 kHz*2pi, omega0=0.177*nu, eta=0.1, detuing=0.95*nu initial=|gg,alpha=sqrt(2)>')
plt.legend()

plt.figure()
xev=np.linspace(-5,5,200)
plt.contourf(xev,xev,wigner(sol.states[-1].ptrace(2),xev,xev),100)
plt.xlabel('x')
plt.ylabel('p')
plt.title('final motional state')
#%%
order=10
hbar=1
size=200
n=np.sqrt(2)
nu=2*np.pi*200*10**3
omega0=nu*0.177
eta=0.1
detuning=nu*0.95
psi=tensor(basis(2,1),basis(2,1),basis(size,2))
phase=0
trange=np.pi/omega0*200
tlist=np.linspace(0,trange,int(100*trange*omega0/np.pi))#np.arange(0,trange,np.pi/omega0/100)
args={'detuning':detuning,'nu':nu,'phase':phase}

sp1=tensor(s_plus,qeye(2),qeye(size))
sm1=tensor(s_minus,qeye(2),qeye(size))
sp2=tensor(qeye(2),s_plus,qeye(size))
sm2=tensor(qeye(2),s_minus,qeye(size))
a=tensor(qeye(2),qeye(2),destroy(size))
H=[]
H.append([0.5*omega0*(sp1+sp2),np.exp(-1j*detuning*tlist)])
H.append([0.5*omega0*(sp1+sp2),np.exp(1j*detuning*tlist)])
H.append([0.5*omega0*(sm1+sm2),np.exp(-1j*detuning*tlist)])
H.append([0.5*omega0*(sm1+sm2),np.exp(1j*detuning*tlist)])
H.append([0.5j*eta*omega0*(sp1+sp2)*a,np.exp(-1j*(nu-detuning)*tlist)])
H.append([0.5j*eta*omega0*(sp1+sp2)*a,np.exp(-1j*(nu+detuning)*tlist)])
H.append([0.5j*eta*omega0*(sp1+sp2)*a.dag(),np.exp(1j*(nu-detuning)*tlist)])
H.append([0.5j*eta*omega0*(sp1+sp2)*a.dag(),np.exp(1j*(nu+detuning)*tlist)])

H.append([-0.5j*eta*omega0*(sm1+sm2)*a,np.exp(-1j*(nu-detuning)*tlist)])
H.append([-0.5j*eta*omega0*(sm1+sm2)*a,np.exp(-1j*(nu+detuning)*tlist)])
H.append([-0.5j*eta*omega0*(sm1+sm2)*a.dag(),np.exp(1j*(nu-detuning)*tlist)])
H.append([-0.5j*eta*omega0*(sm1+sm2)*a.dag(),np.exp(1j*(nu+detuning)*tlist)])

rho=[]
rhoj=[]


sol=mesolve(H,psi,tlist,[],[tensor(sigmaz(),identity(2),identity(size)),tensor(identity(2),sigmaz(),identity(size)),tensor(identity(2),identity(2),num(size))],args=args,options=options)
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
plt.title('nu=200 kHz*2pi, omega0=0.177*nu, eta=0.1, detuing=0.95*nu initial=|gg,2>')
plt.legend()