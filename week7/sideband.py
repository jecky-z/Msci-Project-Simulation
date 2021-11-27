# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 15:41:37 2021

@author: rz3818
"""

#%%
import numpy as np
import scipy as sp
from scipy import integrate
import matplotlib.pyplot as plt
import math
import timeit

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
        raising=raising+np.matmul(a,b)*np.sqrt(i+1+1) # raising|n> -->> sqrt(n+1)|n+1>
    return raising
def lowering(n):
    lowering=np.zeros([n,n])
    for i in range(n-1):
        a=np.zeros([n,n])
        b=np.zeros([n,n])
        a[i,0]=1
        b[0,i+1]=1
        lowering=lowering+np.matmul(a,b)*np.sqrt(i+1+1) # lowering|n> -->> sqrt(n)|n-1>
    return lowering
def unit(n):
    unit=np.zeros([n,n])
    for i in range(n):
        unit[i,i]=1
    return unit
hbar=1

#%%
# Sideband Scan
d_harmonics=10 # number of motional states
miu=1000*2*np.pi
omega0=miu/25
detuning=miu*1
ita=0.1
order=5

def H_sideband(t,n,m):
    if m==n:
        H=0.5*hbar*omega0*(np.exp(-1j*detuning*t)*np.kron(sigma_plus(),unit(d_harmonics))
                            +np.exp(1j*detuning*t)*np.kron(sigma_minus(),unit(d_harmonics)))
    if m>n:
        l=unit(d_harmonics)
        r=unit(d_harmonics)
        for i in range(int(m-n)):
            l=l.dot(lowering(d_harmonics))
            r=r.dot(raising(d_harmonics))
        H=0.5*hbar*omega0*(ita**(m-n))/math.factorial(m-n)*(np.exp(-1j*(detuning+(m-n)*miu)*t)*np.kron(sigma_plus(),l)
                                                +np.exp(1j*(detuning+(m-n)*miu)*t)*np.kron(sigma_minus(),r))
    if m<n:
        l=unit(d_harmonics)
        r=unit(d_harmonics)
        for i in range(int(n-m)):
            l=l.dot(lowering(d_harmonics))
            r=r.dot(raising(d_harmonics))
        H=0.5*hbar*omega0*(ita**(n-m))/math.factorial(n-m)*(np.exp(1j*(detuning-(n-m)*miu)*t)*np.kron(sigma_minus(),l)
                                                +np.exp(-1j*(detuning-(n-m)*miu)*t)*np.kron(sigma_plus(),r))
    return H

def dxdt(t,state):
    dxdt=-1j*np.matmul(H_sideband(t,0,0),state)/hbar
    for i in range(1,1+order):#d_harmonics):
        dxdt=dxdt-1j*np.matmul(H_sideband(t,0,i),state)/hbar-1j*np.matmul(H_sideband(t,i,0),state)/hbar
    return dxdt

n_scan=1000
c=np.zeros([d_harmonics,n_scan]) # collection of excited population at the end of pi pulse
d=np.zeros([n_scan,2,d_harmonics,100]) # collection of full population of all states during pi pulse
start=timeit.default_timer() 
for i in range(n_scan):
    print([i])
    detuning=miu*(-0.5+i/n_scan)*10 # from -5mu to 5mu
    n=2
    psi=np.zeros(d_harmonics)
    psi[n]=1
    psi=np.kron([1,0],psi)+0j # |g,n>
    #omega1=0.25*np.sqrt((omega0*ita*np.sqrt(int(n)+1))**2+(detuning-miu)**2)
    ts=0
    tf=1*np.pi/omega0
    dt=(tf-ts)/100
    t_eval=np.arange(ts,tf,dt)
    #Fs=2*np.pi/dt # sampling frequency
    #f=np.arange(-0.5*len(t_eval),0.5*len(t_eval))*Fs/(len(t_eval)*np.pi)

    sol=sp.integrate.solve_ivp(dxdt,t_span=[ts,tf],y0=psi,t_eval=t_eval,vectorized=True, rtol=1e-7,atol=1e-7)
    y=sol.y.reshape([2,d_harmonics,-1])
    p=y*np.conj(y)
    d[i]=p
    for j in range(d_harmonics):
        c[j,i]=c[j,i]+p[1,j,-1]
end=timeit.default_timer()
print(end-start)
plt.figure()
#for k in range(d_harmonics):
plt.plot((np.arange(n_scan)*10/n_scan-5),np.sum(c,axis=0))
plt.xlabel('detuning / miu')
plt.title('Probability being in |e>')




















