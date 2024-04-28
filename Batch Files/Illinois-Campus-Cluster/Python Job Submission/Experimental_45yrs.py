
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 11 15:06:20 2021

@author: wathe
"""

#%%import libraries
import mip as mp
import numpy as np
from mip import OptimizationStatus
import matplotlib.pyplot as plt
import pandas as pd
import math
import time as time

#%%Define model, parameters, variables and indices
start1=time.time()

m=mp.Model()

x = m.add_var()

#Time frame
t=24

#Alias of t
tp=t

#Number of different O&M alternatives 
n=3

#Number of roads in the network
l=3

#Parameters
#Cost of technique i at time t
#c=[[100000 for i in range(n)] for j in range (t)]
#c=[[3125, 40625, 3125, 32344, 40625, 130000, 218750, 375000, 625000] for j in range(t)]
c=[100000, 50000, 20000]

#IRI deterioration rate
#Building the IRI deterioration matrix to be used for piecewise linear IRI function propagation

#Commented out for a while

prog=pd.read_excel('IRI_Progression_Experimental.xlsx')
d=np.zeros((l,t,t))

for j in range(l):
    d1=list(prog[j+1])
    for i in range (t):
        d[j,i:,i]=d1[0:t-i]
        
RSL=pd.read_excel('Remaining_Service_Life.xlsx')

RSL1=[12,7,3]

#d=3.5

#IRI improvement if technique i is used
#b=[30 for i in range(n)]
b=[100,50,30]

#Unit cost of resurfacing at time t
#kt=[70040 for i in range(t)]
kt=[(70040+abs(np.random.normal()*5000))/2 for i in range(t)]
#Unit cost of milling at time t
#gt=[70040 for i in range(t)]
gt=[(70040+abs(np.random.normal()*5000))/2 for i in range(t)]
#Monetary budget at time t
#mt=[200000 for i in range(t)]
mt=2000000
#A very large number
L=100000
#The upper and lower limits on IRI
IRImax=200
IRImin=60
#Average annual daily traffic
AADTt=[[10000 for i in range (t)] for k in range (l)]
#Average speed
ve=70

#User cost parameters
#Passenger car paramters
kaPC=0.670
kcPC=2.81e-4
dcPC=2.1860e-01
daPC=2.1757e-03
bPC=-1.6931e-01
pPC=3.3753e-04
#Small truck parameters
kaST=7.68e-01
kcST=1.25e-04
dcST=3.0769e-01
daST=0.0108e-03
bST=-7.3026e+01
pST=1.788e-05
#Medium truck parameters
kaMT=9.18e-01
kcMT=1.33e-04
dcMT=9.7418e-01
daMT=9.2993e-03
bMT=+1.3959e-02
pMT=1.0938e+05
#Large truck parameters
kaLT=1.40e+00
kcLT=1.36e-04
dcLT=2.3900e+00
daLT=1.9225e+04
bLT=2.6432e+02
pLT=0.2782e+04

#Traffic split
PC=0.9
ST=0.03
MT=0.03
LT=0.04

#Average parameters
ka=kaPC*PC+kaST*ST+kaMT*MT+kaLT*LT
kc=kcPC*PC+kcST*ST+kcMT*MT+kcLT*LT
dc=dcPC*PC+dcST*ST+dcMT*MT+dcLT*LT
da=daPC*PC+daST*ST+daMT*MT+daLT*LT
bav=bPC*PC+bST*ST+bMT*MT+bLT*LT
p=pPC*PC+pST*ST+pMT*MT+pLT*LT

#Average energy content in a gallon of fuel (MJ)
EC=137.975
#Decision Variables
x = [[[m.add_var(var_type=mp.BINARY) for i in range(t)] for j in range (n) ] for k in range(l)]

y= [[m.add_var(var_type=mp.CONTINUOUS) for i in range (t)] for k in range (l)]

z=[[m.add_var(var_type=mp.CONTINUOUS) for i in range (t)] for k in range (l)]

IRI=[[m.add_var(var_type=mp.CONTINUOUS) for i in range(t)] for k in range (l)]

Et=[[m.add_var(var_type=mp.CONTINUOUS) for i in range(t)] for k in range(l)]

ytp=[[[m.add_var(var_type=mp.BINARY) for i in range(t)] for j in range(tp)] for k in range(l)]

RI=[[m.add_var(var_type=mp.CONTINUOUS) for i in range (t)]for j in range (l)]

x0=[[m.add_var(var_type=mp.BINARY) for i in range (t)]for j in range (l)]

ztp=[[[m.add_var(var_type=mp.CONTINUOUS) for i in range(t)] for j in range(tp)] for k in range(l)]

alpha1=[[m.add_var(var_type=mp.BINARY) for i in range(t)]for j in range (l)]

alpha2=[[m.add_var(var_type=mp.BINARY) for i in range(t)]for j in range (l)]

alpha3=[[m.add_var(var_type=mp.BINARY) for i in range(t)]for j in range (l)]

xn=[[m.add_var(var_type=mp.BINARY) for i in range(t)]for j in range (l)]

de=[[[m.add_var(var_type=mp.BINARY) for i in range(t)]for j in range (l)]for k in range(n)]

start=time.time()

#%%Lagrangian
lambda1=0.5
#Define Constraints
#I
for i in range(t):
    for j in range(l):
        m+=IRI[j][i]<=IRImax
        
#Eqn11 lower limit on IRI
#I
for i in range(t):
    for j in range(l):
        m+=IRI[j][i]>=IRImin
        
#Eqn4 IRI progression equation
#I
for k in range (l):
    m+=IRI[k][0]==60
    for j in range(1,t):
        m+=IRI[k][j-1]-mp.xsum(b[i]*x[k][i][j] for i in range(n))+RI[k][j]==IRI[k][j]
        
#Eqn7 budget constraint
#I
for j in range(t):
    m+=mp.xsum(mp.xsum(c[i]*x[k][i][j] for i in range(n)) for k in range(l))<=mt
    
#Eqn10 auxiliary variable Et[i] calculation
#Energy in MJ/mi
for i in range(t):
    for k in range(l):
        m+=Et[k][i]==((p/ve+(ka*IRI[k][i]+da)+bav*ve+(kc*IRI[k][i]+dc)*ve**2)*AADTt[k][i]*365)/1000

#Eqn12s
#I
for i in range(1,t):
    for j in range (l):
        m+=RI[j][i]==mp.xsum(ytp[j][i][ii]*d[j][i][ii] for ii in range(tp))#+x0[j][i]*d1[i-1]

#UIUN
# for j in range(l):
#     for i in range(t):
#         for ii in range(i,tp):
#             m+=xn[j][i]>=ytp[j][ii][i]
#I
for j in range (l):
    for i in range(t):
        m+=mp.xsum(ytp[j][i][ii] for ii in range(i+1))==1
               
#UI
for j in range (l):
    for i in range(t):
        m+=xn[j][i]==mp.xsum(x[j][k][i] for k in range(n))

#UIUN
# for j in range (l):
#     for i in range(t):
#         m+=xn[j][i]<=1
#I
for j in range(l):        
    for i in range (t):
        for ii in range(i,tp-1):
            m+=ytp[j][ii+1][i]<=ytp[j][ii][i]
#UI            
#questionable
for j in range (l):
    for i in range(t):
            m+=ytp[j][i][i]==xn[j][i]
            
#Remaining service life deduction constraint
for j in range(l):
    for i in range(t):
        for k in range(n):
            m+=de[k][j][i]<=RSL[j][i]*x[j][k][i]



#%%Define Objective Function
# m.objective= mp.xsum(mp.xsum(c[i][j]*x[k][j][i] for i in range(t)) for j in range(n))+mp.xsum(mp.xsum((kt[i]*y[k][i] for i in range(t)))for k in range(l))+mp.xsum(mp.xsum((gt[i]*z[k][i] for i in range(t)))for k in range (l))+mp.xsum(mp.xsum((Et[k][i]/EC for i in range(t))) for k in range(l))*2.5-mp.xsum(mp.xsum(alpha3[j][i] for i in range(t))for j in range(l))

m.objective= mp.xsum(mp.xsum(c[j]*x[k][j][i] for i in range(t)) for j in range(n))+mp.xsum(mp.xsum((Et[k][i]/EC for i in range(t))) for k in range(l))*2.5-mp.xsum(mp.xsum(mp.xsum(((t-i)/RSL1[k]*de[k][j][i])for k in range(l)) for j in range(n))for i in range(t))
#objfun()
#%%Solve the model
#Model's relative gap
m.max_gap=0.001
#Solve the model
m.optimize()
#Show the solution
status = m.optimize(max_seconds=300)

end1=time.time()
totaltime1=end1-start1

t1=[1,2,3,4,5]

x1=np.zeros(l*n*t)

y1=np.zeros(t)

z1=np.zeros(t)

IRI1=np.zeros(t)

Et1=np.zeros(t)

Mill=np.zeros(t)

OL=np.zeros(t)

ytp1=np.zeros(tp*t)

x01=np.zeros(l*t)

RI1=np.zeros(t)

ztp1=np.zeros(tp*t)

alpha11=np.zeros(t)

alpha22=np.zeros(t)

alpha33=np.zeros(t)

xn1=np.zeros((t,l))


end=time.time()

if status == OptimizationStatus.OPTIMAL:
    print('optimal solution cost {} found'.format(m.objective_value))
elif status == OptimizationStatus.FEASIBLE:
    print('sol.cost {} found, best possible: {} '.format(m.objective_value, m.objective_bound))
elif status == OptimizationStatus.NO_SOLUTION_FOUND:
    print('no feasible solution found, lower bound is: {} '.format(m.objective_bound))
if status == OptimizationStatus.OPTIMAL or status == OptimizationStatus.FEASIBLE:
    print('solution:')
    
    
counter=0
for v in m.vars:
#    if abs(v.x) > 1e-6: # only printing non-zeros
    print('{} : {} '.format(v.name, v.x))
    if counter>=0 and counter<(t*n):
        x1[counter]=v.x
    elif counter>=(l*t*n) and counter<(l*t*n+t):
        y1[counter-l*t*n-1]=v.x
    elif counter>=(l*n*t+l*t) and counter<(l*n*t+l*t+t):
        z1[counter-(l*n*t+l*t)-1]=v.x
    elif counter>=(l*n*t+2*l*t) and counter<=(l*n*t+2*l*t+t):
        IRI1[counter-(l*n*t+2*l*t)-1]=v.x
    elif counter>=(l*n*t+3*l*t) and counter<=(l*n*t+3*l*t+t):
        Et1[counter-(l*n*t+3*l*t)-1]=v.x
    elif counter>=(l*n*t+4*l*t) and counter<=(l*n*t+4*l*t+t*tp):
        ytp1[counter-(l*n*t+4*l*t)-1]=v.x
    elif counter>=(l*n*t+4*l*t+t*tp) and counter<=(l*n*t+4*l*t+t*tp+t):
        RI1[counter-(t*n*l+4*t*l+t*tp)-1]=v.x
    elif counter>=(t*n*l+5*t*l+t*tp*l) and counter<=(t*n*l+5*t*l+t*tp*l+t):
        x01[counter-(t*n*l+5*t*l+t*tp*l)-1]=v.x
    elif counter>=(t*n*l+6*t*l+t*tp*l) and counter<=(t*n*l+6*t*l+t*tp*l+t):
        ztp1[counter-(t*n*l+6*t*l+t*tp*l)-1]=v.x
    elif counter>=(t*n*l+6*t*l+t*tp*l+t*tp) and counter<=(t*n*l+6*t*l+t*tp*l+t*tp+t):
        alpha11[counter-(t*n*l+6*t*l+t*tp*l+t*tp)-1]=v.x
    elif counter>=(t*n*l+6*t*l+t*tp*l+t*tp+t) and counter<=(t*n*l+6*t*l+t*tp*l+t*tp+2*t):
        alpha22[counter-(t*n*l+6*t*l+t*tp*l+t*tp+t)-1]=v.x  
    elif counter>=(t*n*l+6*t*l+t*tp*l+t*tp+2*t) and counter<=(t*n*l+6*t*l+t*tp*l+t*tp+3*t):
        alpha33[counter-(t*n*l+6*t*l+t*tp*l+t*tp+2*t)-1]=v.x
    elif counter>=4545 and counter<=4590:
        xn1[counter-4546]=v.x   
    # elif counter>=(l*n*t+4*l*t+t*t) and counter<=(l*n*t+4*l*t+t*t+t):
    #     RI1[counter-(l*n*t+4*l*t)-1]=v.x
    # elif counter>=(l*n*t+5*l*t+t*t) and counter<=(l*n*t+5*l*t+t*t+t):
    #     x0[counter-(l*n*t+5*l*t)-1]=v.x
    counter=counter+1
IRI1[0]=60


end1=time.time()
totaltime1=end1-start1
#%%Sawtooth Plot for the first link

lifetime=np.zeros(t*3-4)
IRI2=np.zeros(t*3-4)
IRI2[0]=60
for i in range(1,len(lifetime)):
    if i%3==1:
        lifetime[i]=math.ceil(i/3)-0.0001
        IRI2[i]=IRI1[math.ceil(i/3)]
    elif i%3==2:
        lifetime[i]=math.ceil(i/3)
        IRI2[i]=IRI1[math.ceil(i/3)]
    elif i%3==0:
        lifetime[i]=math.ceil(i/3)+0.0001
        IRI2[i]=IRI1[math.floor((i)/3)]
plt.plot(lifetime,IRI2)
plt.xlabel('Time (years)')
plt.ylabel('IRI (in/mi)')
plt.savefig('IRIProgrubber.png',dpi=1200)
plt.show()