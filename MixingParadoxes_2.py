#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 15:48:39 2020
Experiments in paradoxical mixing
Experiment 2 : Do perfect gases mix? 
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

#Set initial distribution by patching gas together
ga1=IG.GP(N=800,L=1,beta=0.1); ga1.steps(M=1000,h=1e-2); ga1.q[:,0]+=-1; ga1.q[:,1]+=1
ga2=IG.GP(N=ga1.N,L=1,beta=0.1); ga2.steps(M=1000,h=1e-2);ga2.q[:,0]+=-1; ga2.q[:,1]+=-1

gb1=IG.GP(N=ga1.N,L=1,beta=4); gb1.steps(M=1000,h=1e-2);gb1.q[:,0]+=1; gb1.q[:,1]+=1
gb2=IG.GP(N=ga1.N,L=1,beta=4) ; gb2.steps(M=1000,h=1e-2);gb2.q[:,0]+=1; gb2.q[:,1]+=-1

#COLD GAS
gb=IG.GP(N=ga1.N*2,L=2,beta=gb1.beta)
gb.p[:ga1.N,:]=gb1.p
gb.p[ga1.N:2*ga1.N,:]=gb2.p
gb.q[:ga1.N,:]=gb1.q
gb.q[ga1.N:2*ga1.N,:]=gb2.q

#HOT GAS
ga=IG.GP(N=ga1.N*2,L=2,beta=ga1.beta)
ga.p[:ga1.N,:]=ga1.p
ga.p[ga1.N:2*ga1.N,:]=ga2.p
ga.q[:ga1.N,:]=ga1.q
ga.q[ga1.N:2*ga1.N,:]=ga2.q

#FULL GAS
g=IG.GP(N=ga.N*2,L=2,beta=2/(1/gb1.beta+1/ga1.beta))
g.p[:ga.N,:]=ga.p
g.p[ga.N:2*ga.N,:]=gb.p
g.q[:ga.N,:]=ga.q
g.q[ga.N:2*ga.N,:]=gb.q

output=os.path.join(FigFolder,'Exp2_N%d' %(ga1.N,))
if SAVEFIG and not os.path.exists(output): 
    os.mkdir(output)
    os.mkdir(os.path.join(output,'FULL'))
    os.mkdir(os.path.join(output,'HC'))
    os.mkdir(os.path.join(output,'HC_mix'))


#%% Define Frame 
fig,axs=newfig(1,3,num='Paradoxical Mixing 2',clear=True,tight_layout=True)
bins_P=np.linspace(-4,4,25)
x=np.linspace(-6,6,1000)
dx=x[2]-x[1]

ax=axs[0] #Positions
ax.set_xlim(-3,3)
ax.set_ylim(-3,3)
ax.set_ylabel('$q_y$')
ax.set_xlabel('$q_x$')
lines0=[]
tmp,=ax.plot([],[],linestyle='solid',linewidth=3,color='k');lines0.append(tmp)
tmp, = ax.plot([],[],'o',alpha=0.3,markersize=5, color='forestgreen',markeredgecolor='forestgreen'); lines0.append(tmp)
tmp, = ax.plot([],[],'o',alpha=1,markersize=7, color='k',markeredgecolor='k'); lines0.append(tmp)
tmp, = ax.plot([],[],'o',alpha=0.3,markersize=5, color='firebrick',markeredgecolor='firebrick'); lines0.append(tmp)
tmp, = ax.plot([],[],'o',alpha=0.3,markersize=5, color='mediumblue',markeredgecolor='mediumblue'); lines0.append(tmp)

ax=axs[1] # Energy
ax.set_xlim(0,300)
ax.set_ylim(0,2)
ax.axvline(0,color='firebrick',linestyle='dotted',linewidth=2)
ax.set_ylabel('$E / (N k_BT)$')
ax.set_xlabel('$t$')
lines1=[]
tmp,=ax.plot([],[],'o',markersize=10,markerfacecolor='forestgreen',markeredgecolor='k',label='Kinetic');lines1.append(tmp)
tmp,=ax.plot([],[],'o',markersize=10,markerfacecolor='k',markeredgecolor='k',label='Walls');lines1.append(tmp)
tmp,=ax.plot([],[],'o',markersize=10,markerfacecolor='firebrick',markeredgecolor='k');lines1.append(tmp)
tmp,=ax.plot([],[],'o',markersize=10,markerfacecolor='mediumblue',markeredgecolor='k');lines1.append(tmp)

ax.legend()

def reinit_plot(show_qp=False):
    for lines in lines0+lines1: lines.set_data([],[])
    if show_qp: axs[0].set_ylabel('$p_x$')
    return lines0+lines1

def update_plot(t,show_pdf=True,M=50,h=1e-2,show_hc=True,show_qp=False):
    if g.t>=0:g.L=2        
    #PLOT POSITIONS
    if show_qp:
        if show_hc:
            lines0[3].set_data(g.q[:g.N//2,0],g.p[:g.N//2,0])
            lines0[4].set_data(g.q[g.N//2:,0],g.p[g.N//2:,0])
        else:
            lines0[1].set_data(g.q[:,0],g.p[:,0])
            lines0[2].set_data(g.q[1,0],g.p[1,0])
    else:        
        lines0[0].set_data([-g.L,g.L,g.L,-g.L,-g.L],[-g.L,-g.L,g.L,g.L,-g.L])
        if show_hc:
            lines0[3].set_data(g.q[:g.N//2,0],g.q[:g.N//2,1])
            lines0[4].set_data(g.q[g.N//2:,0],g.q[g.N//2:,1])
        else:
            lines0[1].set_data(g.q[:,0],g.q[:,1])
            lines0[2].set_data(g.q[1,0],g.q[1,1])

    #PLOT ENERGY        
    axs[1].plot(g.t, g.Ekin*g.beta/g.N,'o',color='forestgreen',markersize=1)
    axs[1].plot(g.t, g.Epot*g.beta/g.N,'ok',markersize=1)
    lines1[0].set_data(g.t,g.Ekin*g.beta/g.N)
    lines1[1].set_data(g.t,g.Epot*g.beta/g.N)
    if show_hc:
        axs[1].plot(ga.t, ga.Ekin*g.beta/ga.N,'o',color='firebrick',markersize=1)
        axs[1].plot(gb.t, gb.Epot*g.beta/g.N,'o',color='mediumblue',markersize=1)
        lines1[2].set_data(ga.t, ga.Ekin*g.beta/ga.N)
        lines1[3].set_data(gb.t, gb.Ekin*g.beta/gb.N)

    #PLOT PDF        
    if show_pdf or t==0:
        ax=axs[2]
        ax.clear()
        ax.set_xlim(-3,3)
        ax.set_ylim(0,1)
        if show_hc:
            ax.plot(x,ga.theo(x),linewidth=2,linestyle='dashed',color='firebrick',label='Gibbs  $T= %0.2f$' %(1/ga.beta))
            ax.plot(x,gb.theo(x),linewidth=2,linestyle='dashed',color='mediumblue',label='Gibbs  $T= %0.2f$' %(1/gb.beta))
            tmpa=ga.p.flatten(); tmpa.shape=2*ga.N,1
            tmpb=gb.p.flatten(); tmpb.shape=2*gb.N,1
            axs[2].hist(np.concatenate((tmpa,tmpb),axis=1),bins=bins_P,alpha=0.5,density=True,label=['hot','cold'],color=['firebrick','mediumblue'])
        else:
            ax.plot(x,g.theo(x),linewidth=2,linestyle='dashed',color='k',label='Gibbs  $T= %0.2f$' %(1/g.beta))
            axs[2].hist(g.p,bins=bins_P,alpha=0.5,density=True,label=['$p_x$','$p_y$'])
        ax.legend()
       
    g.steps(M=M,h=h);g.update()
    ga.steps(M=M,h=h);gb.update()
    gb.steps(M=M,h=h);gb.update()

    return lines0+lines1

#%% Set-up for the experiment: We run each gas a,b, full until TMAX
TMAX=100
t=0
pa,qa,pb,qp,p,q=ga.p.copy(),ga.q.copy(),gb.p.copy(),gb.q.copy(),g.p.copy(),g.q.copy()
axs[1].set_xlim(0,TMAX)

while g.t<TMAX:
    update_plot(t,show_hc=False)
    fig.canvas.draw()
    fig.canvas.flush_events()
    if SAVEFIG:
        name=os.path.join(output,'FULL','im%04d.png' %(t,))
        fig.savefig(name)
    t+=1

reinit_plot()
ga.p,ga.q,gb.p,gb.q,g.p,g.q=pa,qa,pb,qp,p,q
ga.t=0
gb.t=0
g.t=0
t=0

while g.t<TMAX:
    update_plot(t,show_hc=True)
    fig.canvas.draw()
    fig.canvas.flush_events()
    if SAVEFIG:
        name=os.path.join(output,'HC','im%04d.png' %(t,))
        fig.savefig(name)
    t+=1
        
reinit_plot(show_qp=True)
axs[0].axvline(-g.L,linewidth=2,color='k')
axs[0].axvline(g.L,linewidth=2,color='k')
axs[0].set_ylim(-6,6)

ga.p,ga.q,gb.p,gb.q,g.p,g.q=pa,qa,pb,qp,p,q
ga.t=0
gb.t=0
g.t=0
t=0
while g.t<TMAX:
    update_plot(t,show_hc=True,show_qp=True)
    fig.canvas.draw()
    fig.canvas.flush_events()
    if SAVEFIG:
        name=os.path.join(output,'HC_mix','im%04d.png' %(t,))
        fig.savefig(name)
    t+=1
