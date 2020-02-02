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

if not os.path.exists('Ehrenfest'): os.mkdir('Ehrenfest')

%matplotlib qt
#%% Entropy of the Ehrenfest model
#The case of distingushable particles
fig,ax=newfig(1,1,num='SB',clear=True,tight_layout=True)
SB = lambda  x : -(x*np.log(x)+(1-x)*np.log(1-x))

SB2= lambda x,N : np.log(scp.special.gamma(N+1))-\
            np.log(scp.special.gamma(N-x+1))-np.log(scp.special.gamma(x+1))

ax.set_xlabel('$a$')
ax.set_ylabel('$S_B/N$')
ax.set_xlim(0,1)
ax.set_ylim(0,1)

ax.axhline(0,linestyle='dotted',color='k')

for N in [10,100]:
    x=np.linspace(0,N,N)
    ax.plot(x/N,SB2(x,N)/N,color='k',label='N=%d' %(N,),linewidth=np.log10(N))

x=np.linspace(0,1,1000)
ax.plot(x,SB(x),color='firebrick',linewidth=5,label='Stirling',alpha=0.6)

ax.legend()

name=os.path.join('Ehrenfest','entropy.pdf')
fig.savefig(name)
#%% Definition of the Ehrenfest model
class Ehrenfest:
    def __init__(self,N=10):
        self.N=N #Number of states
        self.state=1+np.zeros(self.N) #1 : white -1: black
        self.labels=np.arange(self.N)
        self.A= self.compute_A()
        self.B=self.N-self.A
        self.x=(np.random.rand(N)-0.5)*1.9
        self.y=(np.random.rand(N)-0.5)*1.9
        
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
Nit=500#*ehr.N
Time=np.arange(Nit)
A=np.zeros(Nit)

close('all')
fig,axs=newfig(1,2,num='N=%d' %(ehr.N,),clear=True,tight_layout=True)

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
tmp, = ax.plot([],[],'o',color='royalblue',markersize=7,alpha=0.7);lines.append(tmp)
tmp, = ax.plot([],[],'o',color='firebrick',markersize=7,alpha=0.7);lines.append(tmp)

ax.set_xlim(-1.5,1.5)
ax.set_ylim(-1.5,1.5)
ax.plot([-1,1,1,-1,-1],[-1,-1,1,1,-1],linestyle='solid',linewidth=3,color='k',alpha=0.5)


ax.set_xlabel('x')
ax.set_ylabel('y')

def update(i):
    if i>0:
        ehr.update()

    A[i]=ehr.A/ehr.N
    lines[0].set_data(i,ehr.A/ehr.N)
    lines[1].set_data(Time[:i],A[:i])

    x=ehr.x[ehr.state>0]
    y=ehr.y[ehr.state>0]
    lines[2].set_data(x,y)    
    x=ehr.x[ehr.state<0]
    y=ehr.y[ehr.state<0]
    lines[3].set_data(x,y)    
    return lines

#%
output=os.path.join('Ehrenfest','%d' %(ehr.N,))
SAVEFIG=False
if not os.path.exists(output):
    os.mkdir(output)
for i in range(Nit) :
    update(i)
    fig.canvas.draw()
    fig.canvas.flush_events()
    if SAVEFIG:
        name=os.path.join(output,'im%04d.png' %(i,))
        fig.savefig(name)
