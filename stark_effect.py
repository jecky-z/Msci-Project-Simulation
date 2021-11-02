# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 21:14:30 2021

@author: jecky
"""

#%%
import numpy as np
import scipy as sp
from scipy import integrate
import matplotlib.pyplot as plt
import math

def sigma_plus():
    sigma_plus=np.array([[0,1],[0,0]])
    return sigma_plus
def sigma_minus():
    sigma_minus=np.array([[0,0],[1,0]])
    return sigma_minus
def sigma_x():
    sigma_x=np.array([[0,1],[1,0]])
    return sigma_x
def raising(n):
    raising=np.zeros([n,n])
    for i in range(n-1):
        a=np.zeros([n,n])
        b=np.zeros([n,n])
        a[i+1,0]=1
        b[0,i]=1
        raising=raising+np.matmul(a,b)*np.sqrt(i+1)
    return raising
def lowering(n):
    lowering=np.zeros([n,n])
    for i in range(n-1):
        a=np.zeros([n,n])
        b=np.zeros([n,n])
        a[i,0]=1
        b[0,i+1]=1
        lowering=lowering+np.matmul(a,b)*np.sqrt(i+1)
    return lowering
def unit(n):
    unit=np.zeros([n,n])
    for i in range(n):
        unit[i,i]=1
    return unit
hbar=1
#%%
omega0=4*np.pi
detuning=1*np.pi
psi=np.array([1,0])+0j
detunings=100*np.pi
omega1=0.5*np.sqrt((omega0)**2+(detuning)**2)
ts=0
tf=np.pi/omega0
dt=(tf-ts)/1000
t_eval=np.arange(ts,tf,dt)


def H_eff(t): # not in Lamb-Dicke Regime
    H=0.5*hbar*omega0*((np.exp(-1j*detuning*t)+np.exp(-1j*detunings*t))*sigma_plus()
                        +(np.exp(1j*detuning*t)+np.exp(1j*detunings*t))*sigma_minus())
    return H

def dxdt(t,state):
    dxdt=-1j*np.matmul(H_eff(t),state)/hbar
    return dxdt
n_scan=1000
c=[]
for i in range(n_scan):
    print(i)
    detuning=(0.5-i/n_scan)*omega0*10

    t_eval=np.arange(ts,tf,dt)
    sol=sp.integrate.solve_ivp(dxdt,t_span=[ts,tf],y0=psi,t_eval=t_eval,vectorized=True, rtol=1e-7,atol=1e-7)
    Fs=2*np.pi/dt # sampling frequency
    f=np.arange(-0.5*len(t_eval),0.5*len(t_eval))*Fs/(len(t_eval)*2*np.pi)
    p=sol.y*np.conj(sol.y)
    c.append(max(p[1])-min(p[1]))
plt.figure()
plt.plot((0.5-np.arange(n_scan)/n_scan)*10,c)
plt.xlabel('detuning/omega')
plt.title('Population of |e>')
#%%
plt.figure()
plt.subplot(1,2,1)
plt.plot(sol.t,p[0],label='|g>')
plt.plot(sol.t,p[1],label='|e>')
yfft=abs(np.fft.fftshift(np.fft.fft(p[1])))
plt.subplot(1,2,2)
plt.plot(f,yfft,label='|e>')
#plt.xlim([-2*omega1,2*omega1])
plt.subplot(1,2,1)
plt.title('Probability vs t')
plt.legend()
plt.subplot(1,2,2)
plt.title('Frourier Plane')
plt.xlabel('Angular Frequency/pi')
plt.legend()























