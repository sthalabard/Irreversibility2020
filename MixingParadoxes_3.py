#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 15:48:39 2020
Experiments in paradoxical mixing
This time non paradoxical
Experiment 3 : Expansion of a colliding gas 
@author: sthalabard
"""

#%% Startup File
startup_file='init_python.py'
exec(open(startup_file).read())
%matplotlib qt

import CollidingGas as CG


SAVEFIG=True
FigFolder='CollidingGas'
if SAVEFIG and not os.path.exists(FigFolder): os.mkdir(FigFolder)

#Set initial distribution
g=CG.gas(N=40,Lx=1,Ly=1,beta=1,epsilon=1e-3,d0=1e-3,perfect=False)
g.t=-10

output=os.path.join(FigFolder,'Exp3_N%d' %(g.N,))
if SAVEFIG and not os.path.exists(output): 
    os.mkdir(output); os.mkdir(os.path.join(output,'FULL')); os.mkdir(os.path.join(output,'REVERSED'))

#%% Define Frame 
fig,axs=newfig(1,1,num='Non Paradoxical Mixing 3',clear=True,tight_layout=True)
bins_P=np.linspace(-4,4,25)
x=np.linspace(-6,6,1000)
dx=x[2]-x[1]

ax=axs #Positions
ax.set_xlim(-2.1,2.1)
ax.set_ylim(-2.1,2.1)
ax.set_ylabel('$q_y$')
ax.set_xlabel('$q_x$')
lines0=[]
ax.plot([-g.Lx,g.Lx,g.Lx,-g.Lx,-g.Lx],[-g.Ly,-g.Ly,g.Ly,g.Ly,-g.Ly],linestyle='dotted',linewidth=3,color='k',alpha=0.5)

tmp,=ax.plot([],[],linestyle='solid',linewidth=3,color='k');lines0.append(tmp)
tmp, = ax.plot([],[],'o',alpha=0.6,markersize=6, color='mediumblue',markeredgecolor='black'); lines0.append(tmp)
tmp, = ax.plot([],[],'o',alpha=1,markersize=12, color='firebrick',markeredgecolor='k'); lines0.append(tmp)
tmp = ax.text(0.05, 0.9, '', transform=ax.transAxes)
lines0.append(tmp)


def update_plot(t,show_pdf=True,M=10,h=5e-3,showE=True):
    if g.t>=0:g.Lx=2        
    #PLOT POSITIONS
    lines0[0].set_data([-g.Lx,g.Lx,g.Lx,-g.Lx,-g.Lx],[-g.Ly,-g.Ly,g.Ly,g.Ly,-g.Ly])
    lines0[1].set_data(g.q[:,0],g.q[:,1])
    lines0[2].set_data(g.q[0:2,0],g.q[0:2,1])

    g.steps(M=M,h=h)
    g.update()
    return lines0

# Set-up for the experiment:
# We run until t=0, then reverse the velocities to create a ``paradoxical run'' where the gas spontaneously contract at time 100


t=0
while g.t<0:
    update_plot(t,showE=True)
    fig.canvas.draw()
    fig.canvas.flush_events()
    t+=1
    if SAVEFIG:
        name=os.path.join(output,'FULL','im%04d.png' %(t,))
        fig.savefig(name)

g.Lx=2
g.Ly=2
while g.t<100:
    update_plot(t,showE=True)
    fig.canvas.draw()
    fig.canvas.flush_events()
    t+=1
    if SAVEFIG:
        name=os.path.join(output,'FULL','im%04d.png' %(t,))
        fig.savefig(name)


#%%
g.p=-g.p
g.t=0
while g.t<100:
    update_plot(t)
    fig.canvas.draw()
    fig.canvas.flush_events()
    t+=1
    if SAVEFIG:
        name=os.path.join(output,'REVERSED','im%04d.png' %(t,))
        fig.savefig(name)
