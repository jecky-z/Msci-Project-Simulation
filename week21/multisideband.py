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

def Multitone(tlist,nu,omega,eta,detuning,phase,size,order=5,envelope=1,a_error=0):
    H=[]
    for i in range(len(detuning)):
        H += H_sideband(tlist,nu,omega[i],eta,detuning[i]+a_error,phase[i],size,[2,0],order,envelope)
        H += H_sideband(tlist,nu,omega[i],eta,-detuning[i]+a_error,phase[i],size,[2,0],order,envelope)
        H += H_sideband(tlist,nu,omega[i],eta,detuning[i]+a_error,phase[i],size,[2,1],order,envelope)
        H += H_sideband(tlist,nu,omega[i],eta,-detuning[i]+a_error,phase[i],size,[2,1],order,envelope)
    return H

def D(k,eta,size,order=10):
    D=0
    for i in range(order):
        a=(1j*eta)**(2*i+k)
        b=create(size)**(i+k)/math.factorial(i+k)
        c=destroy(size)**i/math.factorial(i)
        D += a*b*c
    return D
    
def Multisideband(tlist,nu,omega,eta,detunings,phase,size,order=2,envelope=1,a_error=0):
    sy=tensor(sigmay(),qeye(2))+tensor(qeye(2),sigmay())
    H=[]
    for i in range(len(detunings)):
        for k in range(0,order+1):
            d=D(k,eta,size)*(1j)**k
            d1=D(k,-eta,size)*(1j)**k
            H.append([-1/eta*omega[i]*tensor(sy,d), envelope*np.exp(1j*((k)*nu-detunings[i])*tlist)*np.exp(1j*a_error*tlist)])
            H.append([-1/eta*omega[i]*tensor(sy,d1), envelope*np.exp(1j*((k)*nu+detunings[i])*tlist)*np.exp(1j*a_error*tlist)])
            H.append([-1/eta*omega[i]*tensor(sy,d.dag()), envelope*np.exp(-1j*((k)*nu-detunings[i])*tlist)*np.exp(1j*a_error*tlist)])
            H.append([-1/eta*omega[i]*tensor(sy,d1.dag()), envelope*np.exp(-1j*((k)*nu+detunings[i])*tlist)*np.exp(1j*a_error*tlist)])
    return H

def Multisideband_0(tlist,nu,omega,eta,detunings,phase,size,order=10,envelope=1,a_error=0):
    sy=tensor(sigmay(),qeye(2))+tensor(qeye(2),sigmay())
    H=[]
    for i in range(len(detunings)):
        d=tensor(sy,D(i+1,eta,size,order=50))*(1j)**i
        H.append([-1/eta*omega[i]*d, envelope*np.exp(1j*((i+1)*nu-detunings[i])*tlist)*np.exp(1j*a_error*tlist)])
        H.append([-1/eta*omega[i]*d.dag(),envelope*np.exp(-1j*((i+1)*nu-detunings[i])*tlist)*np.exp(1j*a_error*tlist)])
    return H

def MSMT_0(tlist,nu,omega,eta,detunings,phase,size,sideband=2,tone=1,order=10,envelope=1,a_error=0):
    sy=tensor(sigmay(),qeye(2))+tensor(qeye(2),sigmay())
    H=[]
    for i in range(sideband):
        for j in range(tone):
            d=tensor(sy,D(i+1,eta,size,order=10))*(1j)**(i+1)
            H.append([-1/eta*omega[tone*i+j]*d, envelope*np.exp(1j*((i+1)*nu-detunings[tone*i+j])*tlist)*np.exp(1j*a_error*tlist)])
            H.append([-1/eta*omega[tone*i+j]*d.dag(),envelope*np.exp(-1j*((i+1)*nu-detunings[tone*i+j])*tlist)*np.exp(1j*a_error*tlist)])
        
    return H


#%%
# Multi Sideband 
order=5
atoms=2
hbar=1
size=100
n=0
nu=2*np.pi*461*10**3
eta=0.3
omega0=nu*0.079
omegak=0.5*omega0*eta*np.exp(0.5*eta**2)
A=np.array([1,1])#-0.5,1.5,1,0])#-0.4264,1.2792,-0.4264,1.2792])
omega=omegak*A



target=(tensor(basis(2,1),basis(2,1)) - 1j*tensor(basis(2,0),basis(2,0)))/np.sqrt(2)
target=tensor(target,qeye(size))
t_rho=target*target.dag()


error=0
c = np.polynomial.Polynomial([-1, 8*eta**2 + 8, -24*eta**2])
def constraint(omega1,omega2,eta):
    a=4*(3+2*eta**2)*omega2**2 + 24*(1+eta**2)*omega1**2
    b=72*omega1**4 + 9*omega2**2 + 24*omega1*omega2**3 + 48*omega1**2*omega2**2
    f=np.polynomial.Polynomial([3,-a,eta**2*b])
    return np.sqrt(np.min(f.roots()))
delta=omegak/np.sqrt(np.min(c.roots()))#constraint(A[0],A[1],eta)#

detuning=nu-delta

detunings=nu*np.array([1,2])-delta*(np.array([2,1])+error)

phase=[0,0,0,0,0,0]
psi=tensor(basis(2,1),basis(2,1),Qobj(thermal_dm(size,n)*(np.zeros([size,1])+1)))
#tensor(ms.Sy.eigenstates()[1][0],basis(size,0))
#basis(2,1).proj(),basis(2,1).proj()

tw=2*np.pi/nu*0
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
envelope=np.ones(len(tlist))
#1/(1+np.exp((2*tw-tlist)/(tw))) - 1/(1+np.exp((trange+6*tw-tlist)/(tw))) # driving field turns on in 0.01 of 2pi/delta
#H=Multitone(tlist,nu,omega,eta,detunings,phase,size,order)
#H=Multisideband_0(tlist,nu,omega,eta,detunings,phase,size,order,envelope)
H=Multisideband(tlist,nu,omega,eta,detunings,phase,size,order,envelope)
#H=MSMT_0(tlist,nu,omega,eta,detunings,phase,size,2,2,order,envelope)


sx=tensor(sigmax(),identity(2),identity(size))+tensor(identity(2),sigmax(),identity(size))
sy=tensor(sigmay(),identity(2),identity(size))+tensor(identity(2),sigmay(),identity(size))
x=tensor(identity(2),identity(2),(destroy(size)+create(size))/np.sqrt(2))
p=tensor(identity(2),identity(2),1j*(-destroy(size)+create(size))/np.sqrt(2))

sol=mesolve(H,psi,tlist,[],[sy*x,sy*p],args=args,options=options)
#(rr*rr.dag())*x,(rr*rr.dag())*p
F=(sol.states[gate_time].dag()*t_rho*sol.states[gate_time]).full()[0]

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


plt.figure(figsize=[9,9])
plt.plot(sol.expect[0][:gate_time+1],sol.expect[1][:gate_time+1],label='1st cycle')
plt.plot(sol.expect[0][gate_time:2*gate_time+1],sol.expect[1][gate_time:2*gate_time+1],label='2nd cycle')
plt.xlabel('<Jy*x>')
plt.ylabel('<Jy*p>')
plt.title('nu=%d kHz*2pi, omega0=%s*nu, eta=%s'%(nu/2000/np.pi,omega0/nu,eta))
plt.legend()
print(sol.expect[0][gate_time],sol.expect[1][gate_time])
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


target=(tensor(basis(2,1),basis(2,1)) - 1j*tensor(basis(2,0),basis(2,0)))/np.sqrt(2)
target=tensor(target,qeye(size))
t_rho=target*target.dag()
F=[]
num=0
options=qutip.solver.Options(atol=1e-12,rtol=1e-12,nsteps=1e6,store_states=True)
loop=np.linspace(0.001,0.01,100)
sideband=1
tone=1
rwa=1
s_error=0
a_error=0
save=True
def constraint(omega1,omega2,eta):
    a=4*(3+2*eta**2)*omega2**2 + 24*(1+eta**2)*omega1**2
    b=72*omega1**4 + 9*omega2**2 + 24*omega1*omega2**3 + 48*omega1**2*omega2**2
    f=np.polynomial.Polynomial([3,-a,eta**2*b])
    return np.sqrt(np.min(f.roots()))


for sideband in [1,2]:
    for tone in [1,2]:
        for o in loop:
            omega0=nu*o*0.5*eta*np.exp(0.5*eta**2)
            print(sideband,tone,o)
            if sideband == 1:
                if tone == 1:
                    omega=np.array([1,1])*omega0
                    delta=4*omega0
                    detuning=nu-delta
                    detunings=np.array([nu-(1+s_error)*delta])
                if tone == 2:
                    omega=np.array([-0.144*4,0.288*4])*omega0
                    delta=4*omega0
                    detuning=nu-delta
                    detunings=np.array([1,1])*nu-(np.array([1,2])+s_error)*delta
            if sideband == 2:
                if tone == 1:
                    omega=np.array([1,1])*omega0
                    c = np.polynomial.Polynomial([-1, 8*eta**2 + 8, -24*eta**2])
                    delta=omega0/np.sqrt(np.min(c.roots()))
                    detuning=nu-delta
                    detunings=np.array([nu-(2+s_error)*delta,2*nu-(1+s_error)*delta])
                if tone ==2:
                    omega=omega0*np.array([-0.5,1.5,1,0])#-0.4264,1.2792,-0.4264,1.2792])
                    c = np.polynomial.Polynomial([-1, 8*eta**2 + 8, -6*eta**2])
                    delta=omega0/np.sqrt(np.min(c.roots()))#constraint(omega[0]/omega0,omega[1]/omega0,eta)
                    detuning=nu-delta
                    detunings=nu*np.array([1,1,2,2])-delta*(np.array([2,6,1,1])+s_error)
                
            phase=[0,0,0,0,0]
            psi=tensor(basis(2,1),basis(2,1),thermal_dm(size,n)*Qobj(1+np.zeros([size,1])))
            
            
            
            tw=2*np.pi/nu*0
            trange=np.pi*(2)/delta
            tstep=np.pi/delta/4000
            tlist=np.arange(0,trange+tw,tstep)
            turnon=np.where(tlist<=tw/2+tstep)[0][-1]
            tp=np.where(tlist<=trange+tw/2)[0][-1]
            gate_time=np.where(tlist <= trange+tw/2)[0][-1]
            args={'detuning':detuning,'nu':nu,'phase':phase}
            #envelope=np.concatenate([np.sin(np.pi*tlist[:turnon]/tw)**2,1+0*tlist[turnon:tp],np.sin(np.pi*(tlist[tp:]-trange)/tw)**2])
            envelope=np.ones([len(tlist)])
            
            if rwa == 1:
                H=Multisideband_0(tlist,nu,omega,eta,detunings,phase,size,order,envelope,a_error*delta)
                if tone > 1:
                    H=MSMT_0(tlist,nu,omega,eta,detunings,phase,size,sideband,tone,order,envelope,a_error*delta)
            if rwa == 0:
                H=Multisideband(tlist,nu,omega,eta,detunings,phase,size,order,envelope,a_error*delta)
            sol=mesolve(H,psi,tlist,[],[],args=args,options=options)
            
            F.append((sol.states[gate_time].dag()*t_rho*sol.states[gate_time]).full()[0])
            
            if save:
                qutip.qsave(sol,'%s'%num)
                np.savetxt('%s.txt'%num,[nu/2/np.pi,omega0/nu,eta,detuning/nu,size,n,s_error,sideband,tone,rwa],
                           header='(trap_frequency(/2*pi), omega0(/nu), eta, detuning(/nu), size, initial state,symetric detuning s_error(/delta),sideband,tone,rwa)'
                           )
                num += 1


F=np.reshape(np.array(F),[4,-1])
plt.figure()
plt.plot(loop,F[0],label='1 sideband 1 tone')
plt.plot(loop,F[1],label='1 sideband 2 tone')
plt.plot(loop,F[2],label='2 sideband 1 tone')
plt.plot(loop,F[3],label='2 sideband 2 tone')

plt.xlabel('symmetric detuing error/delta')#Lamb-Dicke parameter')
plt.ylabel('Infidelity')
plt.title('nu=461 kHz*2pi, omega0=%s*nu, eta=%s,n=%s'%(omega0/nu,eta,n))
plt.legend()
if save:
    np.savetxt('fidelity.txt',F)
    np.savetxt('parameters.txt',[loop[0],loop[-1],nu/2/np.pi,omega0/nu,eta,detuning/nu,size,n,s_error,sideband,tone,rwa],
               header='(start,stop,trap_frequency(/2*pi), omega0(/nu), eta, detuning(/nu), size, initial state,symetric detuning s_error(/delta),sideband,tone,rwa)'
               )

#%%
order=5
atoms=2
hbar=1
size=10
n=0
nu=2*np.pi*200*10**3
omega0=nu*0.00017
omega=np.array([1,1])*omega0
eta=0.1

target=(tensor(basis(2,1),basis(2,1)) - 1j*tensor(basis(2,0),basis(2,0)))/np.sqrt(2)
target=tensor(target,qeye(size))
t_rho=target*target.dag()
F=[]
num=0
options=qutip.solver.Options(atol=1e-12,rtol=1e-12,nsteps=1e6,store_states=True)
loop=np.linspace(0.01,0.5,50)
sideband=1
tone=1
rwa=1
s_error=0
a_error=0
save=True
for sideband in [1,2]:
    for tone in [1,2]:
        for a_error in loop:
            
            print(sideband,tone,a_error)

            if sideband == 1:
                if tone == 1:
                    omega=np.array([1,1])*omega0
                    delta=2*eta*omega0
                    detuning=nu-delta
                    detunings=np.array([nu-(1+s_error)*delta])
                if tone == 2:
                    omega=np.array([-0.144*4,0.288*4,-0.144*4,0.288*4])*omega0
                    delta=2*eta*omega0
                    detuning=nu-delta
                    detunings=np.array([1,1])*nu-(np.array([1,2])+s_error)*delta
            if sideband == 2:
                if tone == 1:
                    omega=np.array([1,1])*omega0
                    p = np.polynomial.Polynomial([-1, 8*eta**2 + 8, -24*eta**2])
                    omegak=0.5*omega0*eta*np.exp(0.5*eta**2)
                    delta=omegak/np.sqrt(np.min(p.roots()))
                    detuning=nu-delta
                    detunings=np.array([nu-(2+s_error)*delta,2*nu-(1+s_error)*delta])
                if tone ==2:
                    omega=np.array([-0.4264,1.2792,-0.4264,1.2792])*omega0
                    p = np.polynomial.Polynomial([-1, 8*eta**2 + 8, -24*eta**2])
                    omegak=0.5*omega0*eta*np.exp(0.5*eta**2)
                    delta=omegak/0.352378
                    detuning=nu-delta
                    detunings=nu*np.array([1,1,2,2])-delta*(np.array([2,4,1,3])+s_error)

                
            phase=[0,0,0,0,0]
            psi=tensor(basis(2,1),basis(2,1),thermal_dm(size,n)*Qobj(1+np.zeros([size,1])))
            
            
            
            trange=np.pi*(2)/delta
            tlist=np.linspace(0,trange,int(1000*trange*delta/np.pi))#np.arange(0,trange,np.pi/omega0/100)
            args={'detuning':detuning,'nu':nu,'phase':phase}
            envelope=1#-np.exp(-(tlist/np.pi/2*delta*50)**2) # driving field turns on in 0.02 of 2pi/delta

            H=Multitone(tlist,nu,omega,eta,detunings,phase,size,order,envelope,a_error*delta)
            sol=mesolve(H,psi,tlist,[],[],args=args,options=options)
            
            F.append((sol.states[-1].dag()*t_rho*sol.states[-1]).full()[0])
            



F=np.reshape(np.array(F),[4,-1])
plt.figure()
plt.plot(loop,1-F[0],label='1 sideband 1 tone')
plt.plot(loop,1-F[1],label='1 sideband 2 tone')
plt.plot(loop,1-F[2],label='2 sideband 1 tone')
plt.plot(loop,1-F[3],label='2 sideband 2 tone')

plt.xlabel('Lamb-Dicke parameter')
plt.ylabel('Infidelity')
plt.title('nu=200 kHz*2pi, omega0=%s*nu, eta=%s,n=%s'%(omega0/nu,eta,n))
plt.legend()
if save:
    np.savetxt('fidelity.txt',F)
    np.savetxt('parameters.txt',[loop[0],loop[-1],nu/2/np.pi,omega0/nu,eta,detuning/nu,size,n,s_error,sideband,tone,rwa],
           header='(start,stop,trap_frequency(/2*pi), omega0(/nu), eta, detuning(/nu), size, initial state,symetric detuning s_error(/delta),sideband,tone,rwa)'
           )


#%%
params = {
   'axes.labelsize': 30,
   'font.size': 30,
   'legend.fontsize': 25,
   'xtick.labelsize': 30,
   'ytick.labelsize': 30,
   'figure.figsize': [20,10],
   'lines.linewidth': 1.5
   } 
plt.rcParams.update(params)
F=np.loadtxt('sd_mt/fidelity.txt',dtype=complex)
parameters=np.loadtxt('sd_mt/parameters.txt')
nu=2*np.pi*parameters[2]
omega0=nu*parameters[3]
eta=parameters[4]
detuning=nu*parameters[5]
delta=nu-detuning
size=parameters[6]
n=parameters[7]
s_error=parameters[8]
sideband=parameters[9]
tone=parameters[10]
rwa=parameters[11]
loop=np.linspace(parameters[0],parameters[1],len(F[0]))

fig,ax=plt.subplots(1,2)
ax[0].plot(loop,F[0],label='Standard MS Gate')
ax[0].plot(loop,F[1],label='Symmetric Detuning Correction')
ax[0].plot(loop,F[2],label='Strong Coupling MS Gate')
ax[0].plot(loop,F[3],label='Compound Gate')
ax[0].set(xlabel='symmetric detuing error/$\delta$', ylabel='Fidelity')
ax[0].grid()
ax[0].set_ylim([0.95,1])

F=np.loadtxt('E:\year4/Msci_Project/codes/Week20/eta/high_on2/fidelity.txt',dtype=complex)
parameters=np.loadtxt('eta/parameters.txt')
loop=np.linspace(parameters[0],parameters[1],len(F[0]))
ax[1].plot(loop,1-F[0],label='Standard MS Gate')
ax[1].plot(loop,1-F[1],label='Symmetric Detuning Correction')
ax[1].plot(loop,1-F[2],label='Strong Coupling MS Gate')
ax[1].plot(loop,1-F[3],label='Compound Gate')
ax[1].set(xlabel='Lamb-Dicke Parameter $\eta$')
ax[1].set_yscale('log')
ax[1].set_xscale('log')
ax[1].grid()

#fig.suptitle('Comparison Across Schemes')#Robustness against Symmetric detuning error')
ax[1].legend()
fig.tight_layout(pad=2)
