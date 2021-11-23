# -*- coding: utf-8 -*-
"""
Created on Wed May 26 09:19:50 2021

@author: Eduardo Paluzo-Hidalgo
"""
from scipy.spatial import Delaunay
import numpy as np
import itertools 

def barycentric(delaunay,points,k): # k=id simplex
    n = len(points[0])
    b = delaunay.transform[k,:n].dot(np.transpose(points - delaunay.transform[k,n]))
    return np.c_[np.transpose(b), 1 - b.sum(axis=0)]

def cartesian(delaunay,barycentric,k):
    return barycentric.dot(delaunay.points[delaunay.simplices[k]])


# Auxiliary function: 
# Input: list
# Output: subsets of a given size
def findsubsets(s, n): 
    return np.array(list(itertools.combinations(s, n)))


# Barycentric subdivision of a given simplex.
# Input: A simplex (simplice in spanish) and points location of the vertices (puntos)
# Output: two lists. The first one is composed of indices of the new vertices
# with the simplices to compute it.
# The second list is composed of the locations of the new points.
def barycentric_subd(simplice,puntos):
    i=len(simplice)
    index = len(puntos)-1
    ls = []
    bs = []
    while i!=1:
        s1= findsubsets(simplice,i)#[0],i)
        b1= [np.sum(puntos[s],axis=0)/i for s in s1]
        xs = []
        for b in b1:
            bs.append(b)
            index+=1
            xs.append(index) 
        ls.append(np.array(list(zip(xs,s1))))
        i-=1
    s_ori= findsubsets(simplice,1)
    ls.append(np.array(list(zip(simplice,s_ori))))
    return (np.concatenate(ls),bs)


# Use the output of barycentric subdivision function to provide the maximal
# simplices of the barycentric subdivision.
# Input: output of barycentric subdivision and dimension of the maximal simplices     
def maximal_simp(t,dim):
    t1=t[0]
    lists = findsubsets(t1,dim+1)
    simps = []
    for l in lists:
        if ite_sublist(l[:,1]):
            simps.append(l[:,0])
    return simps
        
#Auxiliary function to check if a list lst1 is contained in other list lst2
def sublist(lst1, lst2):
    return (set(lst1) <= set(lst2))
 
# Generalization of the function of list to a sequence of lists.
def ite_sublist(ts):
    ls = []
    for i in range(len(ts)):
        a = all([sublist(ts[j],ts[i]) for j in range(i+1,len(ts))])
        ls.append(a)
    return all(ls)
        

# Barycentric subdivision of a simplicial complex.
def bar_subd_complex(simp_complex,puntos):
    ss = []
    p = puntos
    dim=len(simp_complex[0])-1
    for i in range(np.shape(simp_complex)[0]):
        b=barycentric_subd(simp_complex[i],p)
        m = maximal_simp(b,dim)
        p = np.concatenate((p,b[1]))
        ss.append(m)
    return((np.array(np.concatenate(ss),dtype=np.int),p))
    



def nn(dataset,p):
    data = dataset[0] # Points
    classes = dataset[1] # Classes
    k = Delaunay(data) # Complex K
    n = len(set(classes))
    c=np.zeros(n)
    d=np.diag(np.ones(n))
    ls = []
    ls.append(c)
    for i in range(n):
        ls.append(d[i])
    l = Delaunay(ls) 
    b=barycentric(k,p,k.find_simplex(p)[0]) 
    image = classes[k.simplices[k.find_simplex(p)]][0] 
    f=b.dot(l.points[image])
    return  barycentric(l,f,0) 



def nn_simplified(dataset,p):
    data = dataset[0] 
    classes = dataset[1] 
    k = Delaunay(data) 
    pairs = zip(classes[k.simplices],k.simplices)
    triangulos=[y for (x,y) in pairs if len(set(x))!=1]
    k.simplices=np.array(triangulos)
    puntos=k.points[np.unique(k.simplices)]
    classes_simp=classes[np.unique(k.simplices)]
    k=Delaunay(puntos)
    n = len(set(classes))
    c=np.zeros(n)
    d=np.diag(np.ones(n))
    ls = []
    ls.append(c)
    for i in range(n):
        ls.append(d[i])
    l = Delaunay(ls) # COMPLEJO L
    b=barycentric(k,p,k.find_simplex(p)[0]) 
    image = classes_simp[k.simplices[k.find_simplex(p)]][0] 
    f=b.dot(l.points[image])
    return  barycentric(l,f,0) 





    

