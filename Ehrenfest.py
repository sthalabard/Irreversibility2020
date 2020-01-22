#!/home/sthalabard/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  5 13:41:25 2020
Entropy increase in Ehrenfest model
@author: sthalabard
"""

#%% Startup File
startup_file='init_python.py'
exec(open(startup_file).read())

%matplotlib qt
#%% Definition of the Ehrenfest model
class Ehrenfest:
    def __init__(self,N=10):
        self.N=N #Number of states
        self.state=1+np.zeros(self.N) #1 : white -1: black
        self.labels=np.arange(self.N)
        self.A= self.compute_A()
        self.B=self.N-self.A
    def update(self):
        i=np.random.randint(low=0,high=self.N)
        self.state[i]=-1*self.state[i]
        self.A= self.compute_A()
        self.B=self.N-self.A
        return None
        
    def compute_A(self):
        A=(self.state.sum() + self.N)*0.5
        return A

#%% Experiment with N balls
ehr=Ehrenfest(N=1000)
Nit=5*ehr.N
Time=np.arange(Nit)
A=np.zeros(Nit)

close('all')
fig,axs=newfig(1,2,num='N=%d' %(ehr.N,),clear=True)

ax=axs[0] #Time Evolution
ax.axhline(0.5,linestyle='--',linewidth=2,color='k')
#ax.plot(0,ehr.A/ehr.N,'o',\
#               color='royalblue',markeredgecolor='w',markersize=15);
ax.set_xlim(-10,Nit)
ax.set_ylim(-0.01,1.01)
ax.set_xlabel('Time')
ax.set_ylabel('fraction of blue balls')

lines=[]
tmp, = ax.plot([],[],'o',color='royalblue',markersize=10); lines.append(tmp)
tmp, = ax.plot([],[],'--',color='black');lines.append(tmp)

ax=axs[1]
tmp, = ax.plot([],[],'o',color='royalblue',markersize=5,alpha=0.5);lines.append(tmp)
tmp, = ax.plot([],[],'o',color='firebrick',markersize=5,alpha=0.5);lines.append(tmp)

ax.set_xlim(0,ehr.N)
ax.set_ylim(-0.5,1.5)

ax.set_xlabel('labels')
ax.set_ylabel('color')

ax.set_yticks([0,1])
ax.set_yticklabels(['B','R'])

def update(i):
    if i>0:
        ehr.update()

    A[i]=ehr.A/ehr.N
    lines[0].set_data(i,ehr.A/ehr.N)
    lines[1].set_data(Time[:i],A[:i])

    x=ehr.labels[ehr.state>0]
    y=ehr.labels[ehr.state>0]*0
    lines[2].set_data(x,y)    
    x=ehr.labels[ehr.state<0]
    y=ehr.labels[ehr.state<0]*0+1
    lines[3].set_data(x,y)    
    return lines
anim = animation.FuncAnimation(fig, update, init_func=None,repeat=False,
                               frames=Nit, interval=10, blit=True)

fig.tight_layout()
fig.show()
