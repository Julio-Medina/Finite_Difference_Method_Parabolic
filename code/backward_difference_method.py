#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 21:19:53 2023

@author: julio
"""
import numpy as np
import pandas as pd
from crout_factorization_generalization import Crout_generalization 

def f(x): # boundary condition function
    return np.sin(np.pi*x)

def backward_difference_method(l,T,     # l is the x endpoint, T is the maximum time t
                               alpha,   # constant in the parabolic partial d.e.
                               N,m,     # integers defining the grid
                               f):      # functions f
    h=(l)/m
    k=(T)/N
    x=np.linspace(0,l,m+1)
    t=np.linspace(0,T,N+1)
    A=np.zeros((m-1,m-1))
    w=np.zeros(m-1)
    Lambda=alpha**2 *k/(h**2)
    w=f(x[1:-1])
    
    for i in range(m-1):# defines Matrix A
        A[i,i]=(1+2*Lambda)
        if i<m-2:
            A[i,i+1]=-Lambda
            A[i+1,i]=-Lambda
            
    for ti in t:# iterates over time to solve Aw^(j+1)=w^{j}      
        print(ti)
        print(w)
        aux=w[:]
        w=Crout_generalization(A,w,m-1)# uses Crout factorization method
    return A,w,x,aux

def u(x): # function for analytical solution, used in error analysis table
    return np.sin(np.pi*x)*np.exp(-np.pi**2 * 0.5)


def error_table(m,x,w,u): # error analysis table
    csv_list=[]
    for i in range(m+1):
        if i==0 or i==m:
            aux=0
        else:
            aux=w[i-1]
        element=[x[i],aux,u(x[i]),abs(aux-u(x[i]))]
        csv_list.append(element)
    column_scheme=['x_i',
                   'w_i,50',
                   'u(x_i,0.5)',
                   '|u(x_i,0.5)-w_i,50|']
    csv_file_df=pd.DataFrame(csv_list, columns=column_scheme)  
    csv_file_df.to_csv('error_table.csv', index=False)
    return csv_file_df


# defines initial value problem with boundary condition
l=1
T=0.5
alpha=1
N=50
m=10    

# aproximates the solution using the backward difference method(see Burden et. al.)
A,w,x,w_aux=backward_difference_method(l,T,alpha,N,m,f)

# error analysis table
error_Table=error_table(m,x,w_aux,u)
LaTeX_table=error_Table.to_latex()
       
#s=np.sin(np.pi*x)*np.exp(-np.pi**2 * 0.5)
    
