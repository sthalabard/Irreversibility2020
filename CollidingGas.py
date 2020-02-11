#!/home/sthalabard/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  5 13:41:25 2020
Experiments in paradoxical mixing
@author: sthalabard
"""

#%% Startup File
startup_file='init_python.py'
exec(open(startup_file).read())

#%% Definition of a 2D gas with collisions
def _dirac(x,epsilon=1e-3):
    #Smoothened Dirac function are here two smoothen the walls!
    return np.sqrt(1/(2*np.pi*epsilon))*np.exp(-x**2/(2*epsilon))

def _dirac_x(x,epsilon=1e-3):
    #Smoothened derivative of the dirac function 
    return -(x/epsilon)*np.sqrt(1/(2*np.pi*epsilon))*np.exp(-x**2/(2*epsilon))

def _r(x,y,epsilon=1e-3):
    out=x**2+y**2
    out[out<=epsilon**2]= 0.5*(out[out<=epsilon**2] + epsilon**2) 
    return np.sqrt(out)



class gas:
    def __init__(self,N=10,Lx=1,Ly=1,beta=1,m=1,d0=0.1,epsilon=1e-3,perfect=False):
        self.N=N #Number of Molecules
        self.Lx=Lx #Box Size
        self.Ly=Ly
        self.m=m
        self.beta=beta
        self.t=0
        self.d0=d0 #Size of the particles
        self.epsilon=epsilon #Regularization
       #initiate the p distribution with gibbs distrbution     
        self.q,self.p=self.generate_gibbs_QP()
        self.perfect=perfect
        self.update()

        print('E/(N k T)=', self.Ekin*beta /N)
        print('T =',1/self.beta)

    def generate_gibbs_QP(self):
        t=0
        Q=np.zeros((self.N,2))
        Q[:,0]=(np.random.rand(self.N)*2 -1)*self.Lx*0.75
        Q[:,1]=(np.random.rand(self.N)*2 -1)*self.Ly*0.75

        P=np.sqrt(self.m/self.beta)*np.random.randn(self.N,2)
        return Q,P

    def theo(self,x):
        return np.sqrt(self.beta/(2*self.m*np.pi))*np.exp(-self.beta*x**2/(2*self.m))

    def update(self):
        self.Ekin,self.Epot =self.compute_H()

    def Uint(self,q=None,t=None): #Define Interaction potential
        if t is None: t=self.t
        if q is None: q=self.q
        out=np.zeros(self.N)
        if not self.perfect:
            for i in range(self.N):
                Ri= _r(q[i,0]-q[:,0], q[i,1]-q[:,1], epsilon=self.epsilon)
                out[i]=_dirac(Ri-self.d0).sum()
        return out.sum()*0.5
    
    def DUint_DQ(self,q=None,t=None): #Define Interaction potential
        if t is None: t=self.t
        if q is None: q=self.q
        out=np.zeros_like(q)

        if not self.perfect:
            for i in range(self.N):
                Ri= _r(q[i,0]-q[:,0], q[i,1]-q[:,1], epsilon=self.epsilon)
                coeff=0*Ri+1
                coeff[Ri<self.epsilon]=0.5
                out[i,0]=(coeff*_dirac_x(Ri-self.d0)*(q[i,0]-q[:,0])/Ri).sum()
                out[i,1]=(coeff*_dirac_x(Ri-self.d0)*(q[i,1]-q[:,1])/Ri).sum()

        return out

    def U(self,q=None, t=None): #Define enclosing potential 
        if t is None: t=self.t
        if q is None: q=self.q
        out=np.zeros_like(q)
        out[:,0]=_dirac(q[:,0]-self.Lx,epsilon=self.epsilon)+(_dirac(q[:,0]+self.Lx,epsilon=self.epsilon))
        out[:,1]=_dirac(q[:,1]-self.Ly,epsilon=self.epsilon)+(_dirac(q[:,1]+self.Ly,epsilon=self.epsilon))

        return out+self.Uint(q,t)


    def DU_DQ(self, q=None, t=None):
        if t is None: t=self.t
        if q is None: q=self.q
        out=np.zeros((self.N,2))
        #Wall
#        out=dirac_x(q-self.L)+dirac_x(q+self.L)
        out=np.zeros_like(q)
        out[:,0]=_dirac_x(q[:,0]-self.Lx,epsilon=self.epsilon)+_dirac_x(q[:,0]+self.Lx,epsilon=self.epsilon)
        out[:,1]=_dirac_x(q[:,1]-self.Ly,epsilon=self.epsilon)+_dirac_x(q[:,1]+self.Ly,epsilon=self.epsilon)

        #Brute force collision
#        for i in range(self.N):
#            r2=(q[i,0]-q[:,0])**2+(q[i,1]-q[:,1])**2+ epsilon
#            r=np.sqrt(r2)
#            out[i,0]+=(((q[i,0]-q[:,0])/r)*dirac_x(r-self.d0)).sum()
#            out[i,1]+=(((q[i,1]-q[:,1])/r)*dirac_x(r-self.d0)).sum()

        return out+self.DUint_DQ(q,t)

    def compute_H(self,t=None,q=None,p=None): #Compute Hamiltonian
        if t is None: t=self.t
        if q is None: q=self.q
        if p is None: p=self.p
        Ekin = p**2/(2*self.m)
        Epot=self.U(q,t)        
        return Ekin.sum(),Epot.sum()

    def DH_DQ(self,t=None, q=None, p=None):
        if p is None: p=self.p
        if q is None: q=self.q
        if t is None: t=self.t
        out=np.zeros((self.N,2))
        out=self.DU_DQ(q,t)
        return out

    def DH_DP(self,t=None, q=None, p=None):
        if p is None: p=self.p
        if q is None: q=self.q
        if t is None: t=self.t
        out=np.zeros((self.N,2))
        out[:]=p[:]/self.m
        return out

    def step(self,h=1e-2,order=1): #Update in time with symplectic step
        if order==1:
            self.p=self.p-h*self.DH_DQ()
            self.q=self.q+h*self.DH_DP()
        else:#Stormer-Verlet
            phalf = self.p -0.5*h* self.DH_DQ()
            self.q= self.q+h*self.DH_DP(p=phalf)
            self.p=phalf-0.5*h*self.DH_DQ(p=phalf)
        self.t+=h

    def steps(self,h=1e-2,order=2,M=100):
        for i in range(M): self.step(h=h,order=order)
