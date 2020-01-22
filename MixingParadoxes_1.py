#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 15:48:39 2020
Experiments in paradoxical mixing
Experiment 1 : Do perfect gases expand? 
@author: sthalabard
"""

#%% Startup File
startup_file='init_python.py'
exec(open(startup_file).read())
%matplotlib qt

import IdealGas as IG

FigFolder='PerfectGas'
if not os.path.exists(FigFolder): os.mkdir(FigFolder)
SAVEFIG=False

#Set initial distribution
g=IG.GP(N=800,L=1,beta=1)
g.t=-50

output=os.path.join(FigFolder,'Exp1_N%d' %(g.N,))
if SAVEFIG and not os.path.exists(output): 
    os.mkdir(output); os.mkdir(os.path.join(output,'FULL')); os.mkdir(os.path.join(output,'REVERSED'))

#%% Define Frame 
fig,axs=newfig(1,3,num='Paradoxical Mixing 1',clear=True,tight_layout=True)
bins_P=np.linspace(-4,4,25)
x=np.linspace(-6,6,1000)
dx=x[2]-x[1]

ax=axs[0] #Positions
ax.set_xlim(-3,3)
ax.set_ylim(-3,3)
ax.set_ylabel('$q_y$')
ax.set_xlabel('$q_x$')
lines0=[]
ax.plot([-g.L,g.L,g.L,-g.L,-g.L],[-g.L,-g.L,g.L,g.L,-g.L],linestyle='dotted',linewidth=3,color='k',alpha=0.5)

tmp,=ax.plot([],[],linestyle='solid',linewidth=3,color='k');lines0.append(tmp)
tmp, = ax.plot([],[],'o',alpha=0.3,markersize=5, color='mediumblue',markeredgecolor='mediumblue'); lines0.append(tmp)
tmp, = ax.plot([],[],'o',alpha=1,markersize=7, color='firebrick',markeredgecolor='k'); lines0.append(tmp)
tmp = ax.text(0.05, 0.9, '', transform=ax.transAxes)
lines0.append(tmp)

ax=axs[1] # Energy
ax.set_xlim(-50,300)
ax.set_ylim(0,2)
ax.axvline(0,color='firebrick',linestyle='dotted',linewidth=2)
ax.set_ylabel('$E / (N k_BT)$')
ax.set_xlabel('$t$')
lines1=[]
tmp,=ax.plot([],[],'o',markersize=10,markerfacecolor='mediumblue',markeredgecolor='k',label='Kinetic');lines1.append(tmp)
tmp,=ax.plot([],[],'o',markersize=10,markerfacecolor='k',markeredgecolor='k',label='Walls');lines1.append(tmp)
ax.legend()

def update_plot(t,show_pdf=False,M=50,h=1e-2):
    if g.t>=0:g.L=2        
    #PLOT POSITIONS
    lines0[0].set_data([-g.L,g.L,g.L,-g.L,-g.L],[-g.L,-g.L,g.L,g.L,-g.L])
    lines0[1].set_data(g.q[:,0],g.q[:,1])
    lines0[2].set_data(g.q[1,0],g.q[1,1])

    #PLOT ENERGY        
    axs[1].plot(g.t, g.Ekin*g.beta/g.N,'o',color='mediumblue',markersize=1)
    axs[1].plot(g.t, g.Epot*g.beta/g.N,'ok',markersize=1)
    lines1[0].set_data(g.t,g.Ekin*g.beta/g.N)
    lines1[1].set_data(g.t,g.Epot*g.beta/g.N)

    #PLOT PDF        
    if show_pdf or t==0:
        ax=axs[2]
        ax.clear()
        ax.set_xlim(-3,3)
        ax.set_ylim(0,1)
        ax.plot(x,g.theo(x),linewidth=2,linestyle='dashed',color='k',label='Gibbs  $T= %0.2f$' %(1/g.beta))
        axs[2].hist(g.p,bins=bins_P,alpha=0.5,density=True,label=['$p_x$','$p_y$'])
        ax.legend()

    g.steps(M=M,h=h)
    g.update()
    return lines0+lines1

#%% Set-up for the experimient:
# We run until t=100, then reverse the velocities to create a ``paradoxical run'' where the gas spontaneously contract at time 100
t=0
while g.t<100:
    update_plot(t)
    fig.canvas.draw()
    fig.canvas.flush_events()
    t+=1
    if SAVEFIG:
        name=os.path.join(output,'FULL','im%04d.png' %(t,))
        fig.savefig(name)
axs[1].axvline(100,color='firebrick',linestyle='dotted',linewidth=2)

#Copy position and impulses at time t=100   
#Reverse velocities to continue the run 
q100,p100=g.q.copy(),g.p.copy()
g.q,g.p=q100.copy(),-p100.copy()
while g.t<300:
    update_plot(t)
    fig.canvas.draw()
    fig.canvas.flush_events()
    t+=1
    if SAVEFIG:
        name=os.path.join(output,'FULL','im%04d.png' %(t,))
        fig.savefig(name)

# Creation of the Paradoxical run     
g.q,g.p=q100.copy(),-p100.copy()
g.t=0
t=0
axs[1].set_xlim(0,100)
while g.t<100:
    update_plot(t)
    fig.canvas.draw()
    fig.canvas.flush_events()
    t+=1
    if SAVEFIG:
        name=os.path.join(output,'REVERSED','im%04d.png' %(t,))
        fig.savefig(name)