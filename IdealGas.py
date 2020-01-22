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

#%% Definition of a 2D perfect gas
def dirac(x,epsilon=0.001):
    #Smoothened Dirac function are here two smoothen the walls!
    return np.sqrt(1/(2*np.pi*epsilon))*np.exp(-x**2/(2*epsilon))

def dirac_x(x,epsilon=0.001):
    return -(x/epsilon)*np.sqrt(1/(2*np.pi*epsilon))*np.exp(-x**2/(2*epsilon))

class GP:
    def __init__(self,N=10,L=1,beta=1,m=1):
        self.N=N #Number of Molecules
        self.L=L #Box Size
        self.m=m
        self.beta=beta
        self.t=0
        
       #initiate the p distribution with gibbs distrbution     
        self.q,self.p=self.generate_gibbs_QP()
        self.update()

        print('E/(N k T)=', self.Ekin*beta /N)
        print('T =',1/self.beta)

    def generate_gibbs_QP(self,dt=1e-3,Teq=5):
        t=0
        Q=np.zeros((self.N,2))
        P=np.sqrt(self.m/self.beta)*np.random.randn(self.N,2)
        return Q,P

    def theo(self,x):
        return np.sqrt(self.beta/(2*self.m*np.pi))*np.exp(-self.beta*x**2/(2*self.m))

    def update(self):
        self.Ekin,self.Epot =self.compute_H()

    def U(self,q=None, t=None): #Define enclosing potential 
        if t is None: t=self.t
        if q is None: q=self.q
 #       return (1/24)*(np.abs(q)/self.L)**24
        return dirac(q-self.L)+dirac(q+self.L)

    def DU_DQ(self, q=None, t=None):
        if t is None: t=self.t
        if q is None: q=self.q
        out=np.zeros((self.N,2))
#        out= (1/self.L)*(q/self.L) **23
        out=dirac_x(q-self.L)+dirac_x(q+self.L)
        return out

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