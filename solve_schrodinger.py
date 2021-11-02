# -*- coding: utf-8 -*-
"""
Created on Sun Oct 10 16:54:07 2021

@author: jecky
"""
#%%
import numpy as np
import scipy as sp
from scipy import integrate
import matplotlib.pyplot as plt

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
# ODE Solver
def dxdt(t,x):
    dxdt=-1j*np.pi*np.exp(-1j*np.pi*t)
    return dxdt

def solve_ode(function,init_value,dt,endtime):
    result=[]
    result.append(init_value)
    t=np.arange(0,endtime,dt)
    for i in range(len(t)-1):
        result.append(result[i]+function(t[i],result[i])*dt)
    return t,np.array(result)

t,x=solve_ode(dxdt,1,0.001,10)
f=np.exp(-1j*np.pi*t)
plt.plot(t,x,label='numerical')
plt.plot(t,f,label='analytical')
sol=sp.integrate.solve_ivp(dxdt,[0,100],[1],t_eval=t)
plt.plot(sol.t,sol.y[0,:],label='scipy')
plt.title('exp(-i*pi*t) vs t')
plt.legend()
plt.figure()
plt.plot(t,x*np.conj(x),label='numerical')
plt.plot(t,f*np.conj(f),label='analytical')
plt.title('probability vs t')
plt.legend()
#plt.figure()
#plt.plot(t,f-x,label='numerical')
#plt.plot(t,f-sol.y.reshape(sol.t.shape),label='scipy')
#plt.title('difference vs t')
#plt.legend()
#%%
# Rabi Oscillation
hbar=1
omega0=10**2*np.pi # in radian per s
detuning=0
psi=np.array([1,0])+0j
ts=0
tf=10
dt=0.001
t_eval=np.arange(ts,tf,dt)
def H_Rabi(t):
    H=-1j*omega0*0.5*(np.exp(-1j*detuning*t)*sigma_plus()+np.exp(1j*detuning*t)*sigma_minus())
    #H_rabi=-1j*omega0*0.5*np.array([[0, np.exp(-1j*detuning*t)],
    #                                    [np.exp(1j*detuning*t), 0]])
    return H
def dxdt(t,state):
    a=np.matmul(H_Rabi(t),state)
    return a

sol=sp.integrate.solve_ivp(dxdt,t_span=[ts,tf],y0=psi,t_eval=t_eval,vectorized=True, rtol=1e-7,atol=1e-7)
p0=sol.y[0,:]*np.conj(sol.y[0,:]) # 
p1=sol.y[1,:]*np.conj(sol.y[1,:])
#%%
Fs=2*np.pi/dt # sampling frequency
f=np.arange(-0.5*len(sol.t),0.5*len(sol.t))*Fs/(len(sol.t)*2*np.pi) # frequency domain
plt.figure()
plt.plot(sol.t,p0,label='ground')
plt.plot(sol.t,p1,label='excited')
plt.legend()
plt.title('Probability vs t')

plt.figure()
y0fft=abs(np.fft.fftshift(np.fft.fft(p0)))
y1fft=abs(np.fft.fftshift(np.fft.fft(p1)))
plt.plot(f,y0fft,label='ground')
plt.plot(f,y1fft,label='excited')
plt.legend()
plt.title('Fourier Plane')

#%%
# Frequency Scan
omega0=100*np.pi # in radian per s
n_scan=100
c=[]
#a0,a1,f0,f1=[],[],[],[]
for i in range(n_scan):
    detuning=10*np.pi*(i-n_scan/2)
    omega1=0.5*np.sqrt(omega0**2+detuning**2)
    ts=0
    tf=1*np.pi/omega1
    dt=(tf-ts)/10000
    t_eval=np.arange(ts,tf,dt)
    Fs=2*np.pi/dt # sampling frequency
    f=np.arange(-0.5*len(t_eval),0.5*len(t_eval))*Fs/(len(t_eval)*np.pi)
    psi=np.array([1,0])
    sol=sp.integrate.solve_ivp(dxdt,t_span=[ts,tf],y0=np.array([1+0j,0+0j]),t_eval=t_eval,vectorized=True, rtol=1e-7,atol=1e-7)
    p0=sol.y[0,:]*np.conj(sol.y[0,:]) # 
    p1=sol.y[1,:]*np.conj(sol.y[1,:])
    y0fft=abs(np.fft.fftshift(np.fft.fft(p0)))
    y1fft=abs(np.fft.fftshift(np.fft.fft(p1)))
    #index0=sp.signal.find_peaks(y0fft,threshold=10)
    #index1=sp.signal.find_peaks(y1fft,threshold=10)
    #f0.append(f[index0[0][int(len(index0[0])/2)+1:]])
    #f1.append(f[index1[0][int(len(index1[0])/2)+1:]])
    #a0.append(y0fft[index0[0][int(len(index0[0])/2)+1:]])
    #a1.append(y1fft[index1[0][int(len(index1[0])/2)+1:]])
    c.append([p0,p1])
    
    plt.scatter(detuning/omega0,p1[-1],c='b')
plt.xlabel('detuning / Omega')
plt.title('Probability being in |e>')
#%%
j=20
p0=c[j][0][:]
p1=c[j][1][:]
ax=plt.subplot(1,2,1)
plt.suptitle('Omega=%f pi, detuning=%f omega'%(omega0/np.pi,(j-n_scan/2)/100))
plt.figure(1)
plt.plot(t_eval,p0,label='ground')
plt.plot(t_eval,p1,label='excited')
plt.legend()
plt.title('Probability vs t')

plt.subplot(1,2,2)
y0fft=abs(np.fft.fftshift(np.fft.fft(p0)))
y1fft=abs(np.fft.fftshift(np.fft.fft(p1)))
plt.plot(f,y0fft,label='ground')
plt.plot(f,y1fft,label='excited')
plt.xlim([-200,200])
plt.legend()
plt.title('Fourier Plane')
plt.xlabel('angular frequency (2 pi rad/s)')
index0=sp.signal.find_peaks(y0fft,threshold=10)
index1=sp.signal.find_peaks(y1fft,threshold=10)
print(f[index0[0][int(len(index0[0])/2)+1:]])
print(f[index1[0][int(len(index1[0])/2)+1:]])

#%%
# Harmonic motion
d_harmonics=4 # number of motional states
omega0=100*np.pi
detuning=0
ita=0.1
psi=np.zeros([2*d_harmonics])+0j
psi[0]=1+0j # |g,0>


def H_carrier(t):
    H=0.5*hbar*omega0*(np.kron(sigma_plus(),unit(d_harmonics))
        +np.kron(sigma_minus(),unit(d_harmonics)))
    return H

def dxdt(t,state):
    dxdt=-1j*np.matmul(H_carrier(t),state)/hbar
    return dxdt

ts=0
tf=1
dt=0.001
t_eval=np.arange(ts,tf,dt)
sol=sp.integrate.solve_ivp(dxdt,t_span=[ts,tf],y0=psi,t_eval=t_eval,vectorized=True, rtol=1e-7,atol=1e-7)
y=sol.y.reshape([2,4,-1])
Fs=2*np.pi/dt # sampling frequency
f=np.arange(-0.5*len(t_eval),0.5*len(t_eval))*Fs/(len(t_eval)*2*np.pi)
yfft=[]
p=y*np.conj(y)
plt.figure()
for i in range(2):
    for j in range(d_harmonics):
        plt.subplot(1,2,1)
        plt.plot(sol.t,p[i,j],label='|%s,%s>'%(i,j))
        a=abs(np.fft.fftshift(np.fft.fft(p[i,j])))
        yfft.append(a)
        plt.subplot(1,2,2)
        plt.plot(f,a,label='|%s,%s>'%(i,j))
plt.legend()
yfft=np.array(yfft).reshape([2,d_harmonics,-1])

#%%
d_harmonics=4 # number of motional states
omega0=100
detuning=10
miu=10
ita=0.01
psi=np.kron([1,0],[1,0,0,0])+0j # |g,n>
n=np.kron(unit(2),raising(4)).dot(np.kron(unit(2),lowering(4)).dot(psi)).dot(psi.T)
omega1=0.25*np.sqrt(((omega0*ita*np.sqrt(int(n)+1)))**2+(detuning-miu)**2)
ts=0
tf=100
dt=(tf-ts)/10000
t_eval=np.arange(ts,tf,dt)


def H_bluesideband(t):
    H=0.5*hbar*omega0*ita*(np.exp(-1j*(detuning-miu)*t)*np.kron(sigma_plus(),lowering(d_harmonics))
                            +np.exp(1j*(detuning-miu)*t)*np.kron(sigma_minus(),raising(d_harmonics)))
    return H

def dxdt(t,state):
    n=np.kron(unit(2),raising(4)).dot(np.kron(unit(2),lowering(4)).dot(psi)).dot(psi.T)
    dxdt=-1j*np.matmul(H_bluesideband(t),state)/hbar
    return dxdt

sol=sp.integrate.solve_ivp(dxdt,t_span=[ts,tf],y0=psi,t_eval=t_eval,vectorized=True, rtol=1e-7,atol=1e-7)
y=sol.y.reshape([2,4,-1])
Fs=2*np.pi/dt # sampling frequency
f=np.arange(-0.5*len(t_eval),0.5*len(t_eval))*Fs/(len(t_eval)*np.pi)
yfft=[]
p=y*np.conj(y)
plt.figure()
for i in range(2):
    for j in range(d_harmonics):
        plt.subplot(1,2,1)
        plt.plot(sol.t,p[i,j],label='|%s,%s>'%(i,j))
        a=abs(np.fft.fftshift(np.fft.fft(p[i,j])))
        yfft.append(a)
        plt.subplot(1,2,2)
        plt.plot(f,a,label='|%s,%s>'%(i,j))
        plt.xlim([-2*omega1,2*omega1])
plt.subplot(1,2,1)
plt.title('Probability vs t')
plt.legend()
plt.subplot(1,2,2)
plt.title('Frourier Plane')
plt.xlabel('Angular Frequency/pi')
plt.legend()
yfft=np.array(yfft).reshape([2,d_harmonics,-1])
#%%
d_harmonics=4 # number of motional states
omega0=100*np.pi
detuning=10*np.pi
miu=10*np.pi
ita=0.01
psi=np.kron([1,0],[0,1,0,0])+0j # |g,n>
n=np.kron(unit(2),raising(4)).dot(np.kron(unit(2),lowering(4)).dot(psi)).dot(psi.T)
omega1=0.25*np.sqrt((omega0*ita*np.sqrt(int(n)))**2+(detuning-miu)**2)/np.pi
ts=0
tf=2/omega1
dt=(tf-ts)/10000
t_eval=np.arange(ts,tf,dt)


def H_redsideband(t):
    H=0.5*hbar*omega0*ita*(np.exp(1j*(detuning-miu)*t)*np.kron(sigma_plus(),raising(d_harmonics))
                        +np.exp(-1j*(detuning-miu)*t)*np.kron(sigma_minus(),lowering(d_harmonics)))
    return H

def dxdt(t,state):
    dxdt=-1j*np.matmul(H_redsideband(t),state)/hbar
    return dxdt

sol=sp.integrate.solve_ivp(dxdt,t_span=[ts,tf],y0=psi,t_eval=t_eval,vectorized=True, rtol=1e-7,atol=1e-7)
y=sol.y.reshape([2,4,-1])
Fs=2*np.pi/dt # sampling frequency
f=np.arange(-0.5*len(t_eval),0.5*len(t_eval))*Fs/(len(t_eval)*np.pi)
yfft=[]
p=y*np.conj(y)
plt.figure()
for i in range(2):
    for j in range(d_harmonics):
        plt.subplot(1,2,1)
        plt.plot(sol.t,p[i,j],label='|%s,%s>'%(i,j))
        a=abs(np.fft.fftshift(np.fft.fft(p[i,j])))
        yfft.append(a)
        plt.subplot(1,2,2)
        plt.plot(f,a,label='|%s,%s>'%(i,j))
plt.subplot(1,2,1)
plt.title('Probability vs t')
plt.legend()
plt.subplot(1,2,2)
plt.title('Frourier Plane')
plt.xlabel('Angular Frequency/pi')
plt.legend()
yfft=np.array(yfft).reshape([2,d_harmonics,-1])

#%%
# Sideband Scan
d_harmonics=4 # number of motional states
omega0=100
detuning=10
miu=10
ita=0.01


def dxdt(t,state):
    dxdt=-1j*np.matmul(H_bluesideband(t),state)/hbar
    return dxdt

n_scan=100
c=[]

for i in range(n_scan):
    detuning=i*n_scan*0.005
    psi=np.kron([1,0],[0,0,1,0])+0j # |g,n>
    n=n=np.kron(unit(2),raising(4)).dot(np.kron(unit(2),lowering(4)).dot(psi)).dot(psi.T)
    omega1=0.25*np.sqrt((omega0*ita*np.sqrt(int(n)+1))**2+(detuning-miu)**2)
    ts=0
    tf=1/omega1
    dt=(tf-ts)/1000
    t_eval=np.arange(ts,tf,dt)
    Fs=2*np.pi/dt # sampling frequency
    f=np.arange(-0.5*len(t_eval),0.5*len(t_eval))*Fs/(len(t_eval)*np.pi)
    sol=sp.integrate.solve_ivp(dxdt,t_span=[ts,tf],y0=psi,t_eval=t_eval,vectorized=True, rtol=1e-7,atol=1e-7)
    y=sol.y.reshape([2,4,-1])
    p=y*np.conj(y)
    
    plt.scatter(detuning/omega0,p[1,1][-1],c='g')
plt.xlabel('detuning / Omega')
plt.title('Probability being in |e>')
















