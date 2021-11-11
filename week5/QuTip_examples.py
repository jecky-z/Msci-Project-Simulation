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
name'trap_frequency(2*pi), omega0(mu), ita, detuning(mu), size, initial state'

#%%
# Rabi Oscillations
hbar=1
omega0=10**3*2*np.pi
detuning=omega0*0
psi=basis(2,0)
s_plus=(sigmax()-1j*sigmay())*0.5
s_minus=(sigmax()+1j*sigmay())*0.5
trange=5*np.pi/omega0
tlist=np.arange(0,trange,trange/1000)
def H_Rabi(t,arg):
    H=0.5*hbar*omega0*(np.exp(-1j*detuning*t)*s_plus+np.exp(1j*detuning*t)*s_minus)
    return H
sol=mesolve(H_Rabi,psi,tlist,[],[sigmaz()],options=options)
plt.figure()
plt.plot(tlist*omega0/np.pi,sol.expect[0])
plt.xlabel('Time (pi/omega)')
plt.title('Expectation Value')

#%%
# Sideband transitions
hbar=1
size=10
n=2
mu=2*np.pi*10**1
omega0=mu/25
ita=0.1
detuning=mu*1
psi=tensor(basis(2,1),basis(size,n))# used (0,1) as ground state, eigenvalue sigma_z of -1
s_plus=(sigmax()+1j*sigmay())*0.5
s_minus=(sigmax()-1j*sigmay())*0.5
trange=100*np.pi/omega0
tlist=np.arange(0,trange,trange/100)
def H_sideband(t,arg=3):
    s=tensor(s_plus,qeye(size)) # 
    H=0.5*hbar*omega0*(np.exp(1j*(detuning)*t)*s+np.exp(-1j*(detuning)*t)*s.dag())
    for i in range(1,arg+1):
        r=tensor(s_plus,destroy(size)**i)
        b=tensor(s_minus,destroy(size)**i)
        c=0.5*hbar*omega0*(ita**i)/math.factorial(i)
        H=H+c*(np.exp(1j*(detuning-i*mu)*t)*b+np.exp(-1j*(detuning-i*mu)*t)*b.dag()
            +np.exp(-1j*(detuning+i*mu)*t)*r+np.exp(1j*(detuning+i*mu)*t)*r.dag())
    return H
sol=mesolve(H_sideband,psi,tlist,[],[tensor(sigmaz(),qeye(size)),tensor(qeye(2),num(size))],args=3,options=options)
plt.figure()
plt.plot(tlist*omega0/np.pi,sol.expect[0])
plt.xlabel('Time (pi/omega)')
plt.title('Expectation Value')
#%%
# Sideband Scan
hbar=1
size=10
n=2
mu=2*np.pi*10**3
omega0=mu/5
ita=0.1
detuning=mu*2
psi=tensor(basis(2,1),basis(size,n)) # used (0,1) as ground state, eigenvalue sigma_z of -1
s_plus=(sigmax()+1j*sigmay())*0.5
s_minus=(sigmax()-1j*sigmay())*0.5
trange=np.pi/omega0*100
tlist=np.arange(0,trange,trange/100)
def H_sideband(t,arg=3):
    s=tensor(s_plus,qeye(size))
    H=0.5*hbar*omega0*(np.exp(1j*(detuning)*t)*s+np.exp(-1j*(detuning)*t)*s.dag())
    for i in range(1,arg+1):
        r=tensor(s_plus,destroy(size)**i)
        b=tensor(s_minus,destroy(size)**i)
        c=0.5*hbar*omega0*(ita**i)/math.factorial(i)
        H=H+c*(np.exp(1j*(detuning-i*mu)*t)*b+np.exp(-1j*(detuning-i*mu)*t)*b.dag()
            +np.exp(-1j*(detuning+i*mu)*t)*r+np.exp(1j*(detuning+i*mu)*t)*r.dag())
    return H
start=timeit.default_timer()
n_scan=100
scan=(np.arange(0,1,1/n_scan)+200)
c=[]
for i in range(n_scan):
    print(i)
    detuning=scan[i]*mu*0.01
    sol=mesolve(H_sideband,psi,tlist,[],[tensor(sigmaz(),qeye(size)),tensor(qeye(2),num(size))],args=5,options=options)
    c.append(sol)

end=timeit.default_timer()
#%%
expect=np.zeros([n_scan,2,100])
states=np.zeros([n_scan,100,20,1])
for i in range(100):
    expect[i]=c[i].expect
    for j in range(100):
        states[i,j]=c[i].states[j].full()*(np.conj(c[i].states[j].full())) # extract individual states
plt.figure()
plt.plot(scan,expect[:,0,-1])
plt.xlabel('detuing/mu')
plt.title('Expectation sigma_z')
plt.figure()
plt.plot(scan,expect[:,1,-1])
plt.xlabel('detuing/mu')
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
# Comparing run time
run_time=[]
for order in [3,5]:
    for size in [4,10]:
        for detuning in [mu*0,mu*1,mu*2,mu*3]:
            print(order,size,detuning)
            hbar=1
            size=10
            d_harmonics=size
            n=2
            mu=2*np.pi*10**3
            omega0=mu/25
            ita=0.1
            psi=tensor(basis(2,1),basis(size,n))
            s_plus=(sigmax()+1j*sigmay())*0.5
            s_minus=(sigmax()-1j*sigmay())*0.5
            trange=np.pi/omega0*1
            tlist=np.arange(0,trange,trange/100)
            start=timeit.default_timer()
            sol=mesolve(H_sideband,psi,tlist,[],[tensor(sigmaz(),qeye(size)),tensor(qeye(2),num(size))],args=order,options=options)
            end=timeit.default_timer()
            run_time.append(end-start)
            
            psi=np.zeros(size)
            psi[n]=1
            psi=np.kron([1,0],psi)+0j # |g,n>
            #omega1=0.25*np.sqrt((omega0*ita*np.sqrt(int(n)+1))**2+(detuning-miu)**2)
            ts=0
            tf=1*np.pi/omega0
            dt=(tf-ts)/100
            t_eval=np.arange(ts,tf,dt)
            start=timeit.default_timer()
            sol=sp.integrate.solve_ivp(dxdt,t_span=[ts,tf],y0=psi,t_eval=t_eval,vectorized=True, rtol=1e-7,atol=1e-7)
            end=timeit.default_timer()
            run_time.append(end-start)

d_harmonics=10 # number of motional states
miu=1000*2*np.pi
omega0=miu/25
detuning=miu*1
ita=0.1
order=5

def H_1(t,n,m):
    if m==n:
        H=0.5*hbar*omega0*(np.exp(-1j*detuning*t)*np.kron(sigma_plus(),unit(d_harmonics))
                            +np.exp(1j*detuning*t)*np.kron(sigma_minus(),unit(d_harmonics)))
    if m>n:
        l=unit(d_harmonics)
        r=unit(d_harmonics)
        for i in range(int(m-n)):
            l=l.dot(lowering(d_harmonics))
            r=r.dot(raising(d_harmonics))
        H=0.5*hbar*omega0*(ita**(m-n))/math.factorial(m-n)*(np.exp(-1j*(detuning-(m-n)*miu)*t)*np.kron(sigma_plus(),l)
                                                +np.exp(1j*(detuning-(m-n)*miu)*t)*np.kron(sigma_minus(),r))
    if m<n:
        l=unit(d_harmonics)
        r=unit(d_harmonics)
        for i in range(int(n-m)):
            l=l.dot(lowering(d_harmonics))
            r=r.dot(raising(d_harmonics))
        H=0.5*hbar*omega0*(ita**(n-m))/math.factorial(n-m)*(np.exp(1j*(detuning+(n-m)*miu)*t)*np.kron(sigma_minus(),l)
                                                +np.exp(-1j*(detuning+(n-m)*miu)*t)*np.kron(sigma_plus(),r))
    return H

def dxdt(t,state):
    dxdt=-1j*np.matmul(H_1(t,0,0),state)/hbar
    for i in range(order):#d_harmonics):
        dxdt=dxdt-1j*np.matmul(H_1(t,0,i),state)/hbar-1j*np.matmul(H_1(t,i,0),state)/hbar
    return dxdt

n_scan=1000
c=np.zeros([d_harmonics,n_scan])
d=np.zeros([n_scan,2,d_harmonics,100])
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


