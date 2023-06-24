#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 29 14:31:32 2022

@author: spri902
"""
# Homework Instruction: Please follow order of this notebook and fill in the codes where commented as TODO
# require numpy module
import numpy as np
from numpy.linalg import norm
from numpy.linalg import solve
import matplotlib.pyplot as plt
# import the supplementary function for homework 1
import sys
sys.path.insert(0, './')
from hw1_supp import *
#Gradient Descent Solver
#Recall the gradient descent algorithm we learned from the class and complete the gradient descent solver.
def optimizeWithGD(x0, func, grad, step_size, tol=1e-6, max_iter=1000, use_line_search=False):
    """
    Optimize with Gradient Descent
    
    input
    -----
    x0 : array_like
        Starting point for the solver.

    func : function
        Input x and return the function value.
    grad : function
        Input x and return the gradient.
    step_size : float or None
        If it is a float number and `use_line_search=False`, it will be used as the step size.
        Otherwise, line search will be used
    tol : float, optional
        Gradient tolerance for terminating the solver.
    max_iter : int, optional
        Maximum number of iteration for terminating the solver.
    use_line_search : bool, optional
        When it is true a line search will be used, other wise `step_size` has to be provided.
        
    output
    ------
    x : array_like
        Final solution
    obj_his : array_like
        Objective function value convergence history
    err_his : array_like
        Norm of gradient convergence history
    exit_flag : int
        0, norm of gradient below `tol`
        1, exceed maximum number of iteration
        2, line search fail
        3, others
    """
    
    # safeguard
    if not use_line_search and step_size is None:
        print('Please specify the step_size or use the line search.')
        return x0, np.array([]), np.array([]), 3
    
    # initial step
    x = np.copy(x0)
    g = grad(x)
    #
    obj = func(x)
    err = norm(g)
    #
    obj_his = np.zeros(max_iter + 1)
    err_his = np.zeros(max_iter + 1)
    #
    obj_his[0] = obj
    err_his[0] = err
    # start iterations
    iter_count = 0
    while err >= tol:
        if use_line_search:
            step_size = lineSearch(x, g, g, func)
        #
        # if line search fail step_size will be None
        if step_size is None:
            print('Gradient descent line search fail.')
            return x, obj_his[:iter_count+1], err_his[:iter_count+1], 2
        #
        # gradient descent step
        #####
        # TODO: with given step_size, complete gradient descent step
        # x = ? 
        #####
        x -= step_size * g
        #
        # update function and gradient
        g = grad(x)
        #
        #
        obj = func(x)
        err = norm(g)
        #
        iter_count += 1
        obj_his[iter_count] = obj
        err_his[iter_count] = err
        #
        # check if exceed maximum number of iteration
        if iter_count >= max_iter:
            print('Gradient descent reach maximum number of iteration.')
            return x, obj_his[:iter_count+1], err_his[:iter_count+1], 1
    #
    return x, obj_his[:iter_count+1], err_his[:iter_count+1], 0

#Newton's Solver
#Recall the Newton's method we learned from the class and complete the Newton's solver.

#Remark: This is a simplified version of Newton's method which do not evolve line search with constant step size 1.
def optimizeWithNT(x0, func, grad, hess, tol=1e-6, max_iter=100):
    """
    Optimize with Newton's Method
    
    input
    -----
    x0 : array_like
        Starting point for the solver.
    func : function
        Input x and return the function value.
    grad : function
        Input x and return the gradient.
    hess : function
        Input x and return the Hessian matrix.
    tol : float, optional
        Gradient tolerance for terminating the solver.
    max_iter : int, optional
        Maximum number of iteration for terminating the solver.
        
    output
    ------
    x : array_like
        Final solution
    obj_his : array_like
        Objective function value convergence history
    err_his : array_like
        Norm of gradient convergence history
    exit_flag : int
        0, norm of gradient below `tol`
        1, exceed maximum number of iteration
        2, others
    """
    # initial step
    x = np.copy(x0)
    g = grad(x)
    H = hess(x)
    #
    obj = func(x)
    err = norm(g)
    #
    obj_his = np.zeros(max_iter + 1)
    err_his = np.zeros(max_iter + 1)
    #
    obj_his[0] = obj
    err_his[0] = err
    # start iteration
    iter_count = 0
    while err >= tol:
        # Newton's step
        #####
        # TODO: complete Newton's step
        # x = ? hint: using solve function
        #####
        x -= solve(H,g)
        #
        # update function, gradient and Hessian
        g = grad(x)
        H = hess(x)
        #
        #
        obj = func(x)
        err = norm(g)
        #
        iter_count += 1
        obj_his[iter_count] = obj
        err_his[iter_count] = err
        #
        # check if exceed maximum number of iteration
        if iter_count >= max_iter:
            print('Gradient descent reach maximum number of iteration.')
            return x, obj_his[:iter_count+1], err_his[:iter_count+1], 1
    #
    return x, obj_his[:iter_count+1], err_his[:iter_count+1], 0

#Simple Test
#Consider a simple test problem:
   # min f(x) = 1/2*||x-b||^2

#We could easily calculate the gradient and Hessian of f



#   gradient(f(x)) = x-b, laplacian(f(x)) = I

#where I is the identity matrix and  we know that x* = b is the solution

#Please complete the solvers above, use the following function as a simple test.

#Remark: Having simple test is always a good idea. It can easily detect the 
#bugs in the code and save a lot of time and emotion energy when use the code 
#in the more complex context.

# create b
b = np.array([1.0, 2.0, 3.0])
# define test function
def test_func(x):
    return 0.5*sum((x - b)**2)
# define test gradient
def test_grad(x):
    return x - b
# define test Hessian
def test_hess(x):
    return np.eye(b.size)

# test gradient descent
x0_gd = np.zeros(b.size)
#
x_gd, obj_his_gd, err_his_gd, exit_flag_gd = optimizeWithGD(x0_gd, test_func, test_grad, 1.0)
# check if solution is correct
err_gd = norm(x_gd - b)
#
if err_gd < 1e-6:
    print('Gradient Descent: OK')
else:
    print('Gradient Descent: Err')
    
# test Newton's method
x0_nt = np.zeros(b.size)
#
x_nt, obj_his_nt, err_his_nt, exit_flag_nt = optimizeWithNT(x0_nt, test_func, test_grad, test_hess)
# check if solution is correct
err_nt = norm(x_nt - b)
#
if err_nt < 1e-6:
    print('Newton: OK')
else:
    print('Newont: Err')


#Logistic Regression
#Recall the logistic regression model from the class,
#            m
#    f(x) = sum(log[1 + exp((<a,x>))] - b<a,x>) + (lambda/2)||x||^2
#           i=1
#Calculate the gradient and Hessian of the function and define 
#them below using the given data.

# fix the random seed to generate same data set
np.random.seed(123)
# set dimension and create the data
m_lgt = 500
n_lgt = 50
A_lgt = 0.3*np.random.randn(m_lgt, n_lgt)
x_lgt = np.random.randn(n_lgt)
b_lgt = sampleLGT(x_lgt, A_lgt)
lam_lgt = 0.1

#Define functions
# define function, gradient and Hessian
def lgt_func(x):
    #####
    # TODO: complete the function
    #####
    return np.sum(np.log(1+np.exp(A_lgt@x))) - b_lgt@A_lgt@x + lam_lgt*x@x/2
#
def lgt_grad(x):
    #####
    # TODO: complete the function
    #####
    return  A_lgt.T@ ((np.exp(A_lgt@x)/(1+np.exp(A_lgt@x))) - b_lgt) + lam_lgt*x
#
def lgt_hess(x):
    #####
    # TODO: complete the function
    #####
    return A_lgt.T @ np.diag( np.exp(A_lgt@x)/(1+np.exp(A_lgt@x))**2 ) @ A_lgt + lam_lgt * np.eye(n_lgt)

# uniform step size for gradient descent or should we use line search?
#####
# TODO: if there is a uniform stepsize, set up step_size_lgt = 
# otherwise set use_line_search=True
#####
step_size_lgt = 1/( np.sqrt(n_lgt)/2*np.linalg.norm(A_lgt,2)**2 + lam_lgt )

#Apply a gradient descent solver
x0_lgt_gd = np.zeros(n_lgt)
x_lgt_gd, obj_his_lgt_gd, err_his_lgt_gd, exit_flag_lgt_gd = \
    optimizeWithGD(x0_lgt_gd, lgt_func, lgt_grad, step_size_lgt)
    
# plot result
fig, ax = plt.subplots(1, 2, figsize=(12,5))
ax[0].plot(obj_his_lgt_gd)
ax[0].set_title('function value')
ax[1].semilogy(err_his_lgt_gd)
ax[1].set_title('norm of gradient')
fig.savefig('img/logistic_gd.pdf')
fig.suptitle('Gradient Descent on Logistic Regression')
plt.show()

#Apply Newton's Solver
x0_lgt_nt = np.zeros(n_lgt)
x_lgt_nt, obj_his_lgt_nt, err_his_lgt_nt, exit_flag_lgt_nt = \
    optimizeWithNT(x0_lgt_nt, lgt_func, lgt_grad, lgt_hess)
    
# plot result
fig, ax = plt.subplots(1, 2, figsize=(12,5))
ax[0].plot(obj_his_lgt_nt)
ax[0].set_title('function value')
ax[1].semilogy(err_his_lgt_nt)
ax[1].set_title('norm of gradient')
fig.savefig('img/logistic_nm.pdf')
fig.suptitle('Newton\'s Method on Logistic Regression')
plt.show()


#Poisson Regression
#Recall the logistic regression model from the class:
#            m
#    f(x) = sum(exp(<a,x>) - b<a,x>) + (lambda/2)||x||^2
#           i=1

#Calculate the gradient and Hessian of the function and define 
#them below using the given data


# fix the random seed to generate same data set
np.random.seed(123)
# set dimension and create the data
m_psn = 500
n_psn = 50
A_psn = 0.3*np.random.randn(m_psn, n_psn)
x_psn = np.random.randn(n_psn)
b_psn = samplePSN(x_psn, A_psn)
lam_psn = 0.1

#Define functions
# define function, gradient and Hessian
def psn_func(x):
    #####
    # TODO: complete the function
    #####
    return np.sum(np.exp(A_psn@x)) - b_psn.T@A_psn@x + lam_psn*x@x/2
#
def psn_grad(x):
    #####
    # TODO: complete the function
    #####
    return A_psn.T@(np.exp(A_psn@x) - b_psn) + lam_psn*x
#
def psn_hess(x):
    #####
    # TODO: complete the function
    #####
    return  A_psn.T*np.exp(A_psn@x)@A_psn +  lam_psn * np.eye(n_psn)
# uniform step size for gradient descent or should we use line search?
#####
# TODO: if there is a uniform stepsize, set up step_size_lgt = 
# otherwise set use_line_search=True
#####

#Apply gradient descent solver
x0_psn_gd = np.zeros(n_psn)
x_psn_gd, obj_his_psn_gd, err_his_psn_gd, exit_flag_psn_gd = \
    optimizeWithGD(x0_psn_gd, psn_func, psn_grad, None, use_line_search=True)
    
# plot result
fig, ax = plt.subplots(1, 2, figsize=(12,5))
ax[0].plot(obj_his_psn_gd)
ax[0].set_title('function value')
ax[1].semilogy(err_his_psn_gd)
ax[1].set_title('norm of gradient')
fig.savefig('img/poisson_gd.pdf')
fig.suptitle('Gradient Descent on Poisson Regression')
plt.show()

#Apply Newton's Solver
x0_psn_nt = np.zeros(n_lgt)
x_psn_nt, obj_his_psn_nt, err_his_psn_nt, exit_flag_psn_nt = \
    optimizeWithNT(x0_psn_nt, psn_func, psn_grad, psn_hess)

# plot result
fig, ax = plt.subplots(1, 2, figsize=(12,5))
ax[0].plot(obj_his_psn_nt)
ax[0].set_title('function value')
ax[1].semilogy(err_his_psn_nt)
ax[1].set_title('norm of gradient')
fig.savefig('img/poisson_nm.pdf')
fig.suptitle('Newton\'s Method on Poisson Regression')
plt.show()
































