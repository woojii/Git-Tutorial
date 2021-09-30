#!/usr/bin/env python
# coding: utf-8

# In[1]:


#import TB model class - atoms in 1d 
from pythtb import *
import matplotlib.pyplot as plt

#specify model
#lattice vector (we only have one)-단위 행렬
lattice=[[1.0]] 
#positons of orbitals (only one)
orb=[[0.0]] 

#dim_k: 상호공간의 차원, 얼마나 많은 방향이 주기적으로 간주되는지
#dim_r: 실제 공간의 차원, 실제 공간 격자 벡터가 몇개 있고 궤도 좌표를 지정하는데 필요한 좌표가 몇개인지
#define the model(dimension k space, dimension real space, lattice vecs, orbital vecs)
mymodel=tb_model(1,1,lattice,orb)

#assign hopping terms(amplitude,iα,jβ,R)
#amplitude- 도약진폭, 실수 또는 복소수 , 단일 숫자로 주어지면 업 및 다운스핀 구성 모두에 대한 도약 진폭
#R - ket 오비탈이 위치한 단위 셀을 가리키는 격자 벡터, dim_r 과 같아야 함
#단단히 바인딩된 오비탈 간의 호핑매개변수
mymodel.set_hop(1.,0,0,[1]) 

#define a path in k-space
# k path to plot in units of reciprocal lattice vecs (-π to π)
path=[[-0.5],[0],[0.5]] 
label=(r'$-\pi$',r'$0$',r'$\pi $')#label k points
numsteps=100 # number of steps between points
# k_path(kpts,nk) kpts:지정된 k 포인트 사이 상호공간에서 경로 보간, nk:플롯만드는데 사용할 총 k point 수
#kpts : full-[0.0,0.5,1.0] 전체 BZ, fullc-[-0.5,0.0,0.5]전체BZ,중앙
kpts=k_path(path,numsteps) 

#solve for eigen energies of hailtonian
#solve_all (k_list) : 주어진 1차원 k벡터 목록에서 긴밀한 결합 모델의 고유값 및 고유벡터를 푼다
evals=mymodel.solve_all(kpts)#solve model

fig=plt.figure()#make a figure object
plt.plot(evals[0])#repeat this line for more bands, up to evals[n]
plt.title('1d chain of atoms')
plt.xlabel('path in k space')
plt.ylabel('band energy')
plt.xticks([0,100,200],label)


# In[34]:


#specifiy model - 2-d SQUARE LATTICE
lattice1=[[1.0,0.0],[0.0,1.0]] # two component lattice vectors
orb1=[[0.0,0.0]]# still one orbital, two components

#define the model
mymodel1=tb_model(2,2,lattice1,orb1)
#x hopping, y hopping
mymodel1.set_hop(-1.,0,0,[1.0,0])
mymodel1.set_hop(-1.,0,0,[0,1.0]) 
     
#define a path in k space to plot (r,x,m,r)
path=[[0.0,0.0],[.5,0],[.5,.5],[0.0,0.0]]
     
label=(r'$-\Gamma $',r'$X$',r'$M $',r'$-\Gamma $')#label k points
numsteps=100 # number of steps between points
kpts=k_path(path,numsteps)

evals=mymodel1.solve_all(kpts)#solve model

fig=plt.figure()#make a figure object
plt.plot(evals[0])#repeat this line for more bands, up to evals[n]
plt.title('2d square lattice of atoms')
plt.xlabel('path in k space')
plt.ylabel('band energy')
plt.xticks([0,100,200,300],label)


# In[38]:


#specify model
#lattice vectors
lat2=[[1.0,0.0],[0.0,1.0]]
#positions of orbitals - 2 orbitals
orb2=[[0.0,0.0],[.5,.5]]

#define the model
my_model=tb_model(2,2,lat2,orb2)

#assign onsite energy for each orbital(onsite_en)
#단단히 결합된 궤도에 대한 현장에너지
my_model.set_onsite([1.0,-1.0])

#assign hopping terms
t=1
t2=1
#x-hopping, y-hopping within sublattice of orbital'0' (A sublattice hopping)
my_model.set_hop(-t,0,0,[1.0,0])
my_model.set_hop(-t,0,0,[0,1.0])
#x-hopping, y-hopping within sublattice of orbital'1' (B sublattice hopping)
my_model.set_hop(-t,1,1,[1.0,0])
my_model.set_hop(-t,1,1,[0,1.0])
#4 inter-sublattice hopping terms, from 0 to 1(last term is R)
my_model.set_hop(-t2,0,1,[0.0,0.0])
my_model.set_hop(-t2,0,1,[-1.0,0.0])
my_model.set_hop(-t2,0,1,[-1.0,-1.0])
my_model.set_hop(-t2,0,1,[0.0,-1.0])

path=[[0.0,0.0],[.5,0],[.5,.5],[0.0,0.0]]
     
label=(r'$-\Gamma $',r'$X$',r'$M $',r'$-\Gamma $')#label k points
numsteps=100 # number of steps between points
kpts=k_path(path,numsteps)

evals=my_model.solve_all(kpts)#solve model

fig=plt.figure()#make a figure object
plt.plot(evals[0])
plt.plot(evals[1])#repeat this line for more bands, up to evals[n]
plt.title('interpenetrating square lattices')
plt.xlabel('path in k space')
plt.ylabel('band energy')
plt.xticks([0,100,200,300],label)

