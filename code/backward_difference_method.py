#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 21:19:53 2023

@author: julio
"""
import numpy as np
import pandas as pd
from crout_factorization_generalization import Crout_generalization 

def f(x):
    return np.sin(3.141592*x)

def backward_difference_method(l,T,     # l is the x endpoint, T is the maximum time t
                               alpha,   # constant in the parabolic partial d.e.
                               N,m,     # integers defining the grid
                               f,u):    # functions f and u 
    h=(l)/m
    k=(T)/N
    x=np.linspace(0,l,m+1)
    t=np.linspace(0,T,N+1)
    A=np.zeros((m-1,m-1))
    w=np.zeros(m-1)
    Lambda=alpha**2 *k/(h**2)
    w=f(x[1:-1])
    
    for i in range(m-1):
        A[i,i]=(1+2*Lambda)
        if i<m-2:
            A[i,i+1]=-Lambda
            A[i+1,i]=-Lambda
    for ti in t:
         
        print(ti)
        print(w)
        w=Crout_generalization(A,w,m-1)
    return A,w,x

A,w,x=backward_difference_method(1,0.5,     # l is the x endpoint, T is the maximum time t
                               1,   # constant in the parabolic partial d.e.
                               50,10,     # integers defining the grid
                               f,f)
        
    
    
s=np.sin(np.pi*x)*np.exp(-np.pi**2 * 0.5)
    
