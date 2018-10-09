###################################################
##                                               ##
## This file is part of the KinBot code v2.0     ##
##                                               ##
## The contents are covered by the terms of the  ##
## BSD 3-clause license included in the LICENSE  ##
## file, found at the root.                      ##
##                                               ##
## Copyright 2018 National Technology &          ##
## Engineering Solutions of Sandia, LLC (NTESS). ##
## Under the terms of Contract DE-NA0003525 with ##
## NTESS, the U.S. Government retains certain    ##
## rights to this software.                      ##
##                                               ##
## Authors:                                      ##
##   Judit Zador                                 ##
##   Ruben Van de Vijver                         ##
##                                               ##
###################################################
import os,sys
import math
import logging
import numpy as np


class BFGS:
    """
    Class to execute a BFGS optimization to a local minimum
    """
    def __init__(self,step_tol = 1E-7, grad_tol = 1E-7, line_tol = 1E-10, inhess = [],
                 max_step = 100, max_lin_step = 1000, use_grad_tol = 1, use_step_tol = 1):
        """
        Initialize the BFGS algorithm
        
        step_tol: tolerance on the minimum step size to continue the optimization (absolute L2 norm of the step)
        
        grad_tol: tolerance on the minimum gradient to continue the optimization (absolute L2 norm of the gradient)
        
        line_tol: tolerance on the alpha in the line search
        
        inhess: initial guess for the hessian (default is identity matrix)
        
        max_step: maximum number of Gradient iterations
        
        max_lin_step: maximum number of iterations in the linear search
        
        use_grad_tol: use the gradient tolerance convergance criterion
        
        use_step_tol: use the step size tolerance convergance criterion
        """
        self.step_tol = step_tol
        self.grad_tol = grad_tol
        self.line_tol = line_tol
        self.inhess = inhess
        self.max_step = max_step
        self.max_lin_step = max_lin_step
        self.use_grad_tol = use_grad_tol
        self.use_step_tol = use_step_tol
        
        if not use_grad_tol and not use_step_tol:
            sys.exit('Cannot execute an optimization if neither the step nor the gradient tolerance can be used')
            
    
    def converged(self,step,grad):
        if np.linalg.norm(step) < self.step_tol and self.use_step_tol:
            return 1
        if np.linalg.norm(grad) < self.grad_tol and self.use_grad_tol:
            return 1
        return 0

    def optimize(self,f,x):
        """
        Optimize a function to the closest local minimum using the BFGS method.
        
        f: a function of x which needs to have
            f.eval(x): function evaluation of f
            f.gradient(x): gradient df/dx in x
        
        x: initial guess
        """
        
        x_i = [x] #intermediate coordinates
        
        if len(self.inhess) == 0:
            H = np.eye(len(x))
        else:
            H = self.inhess

        g = f.gradient(x)
        g_i = [g] #intermediate forces
        p = np.dot(np.linalg.inv(H),-g)
        
        it = 0
        while it < self.max_step:
            ak = self.line_search(f,x,p)
            sk = ak * p
            xk = x + sk
            gk = f.gradient(xk)
            yk = gk - g
            
            if np.dot(sk,yk):
                t1 = np.outer(yk, yk) / np.dot(sk,yk)
            else:
                t1 = np.outer(yk, yk)
            t2 = np.dot(H,sk)
            t3 = np.outer(t2,t2) / np.dot(sk,t2)
            Hk = H + t1 - t3
            
            x_i.append(xk)
            g_i.append(gk)
            
            if self.converged(xk - x, gk):
                return xk, x_i, g_i
            g = gk
            x = xk
            p = np.dot(np.linalg.inv(Hk),-gk)
            H = Hk
            
            it += 1
        
        return xk, x_i, g_i
    
    def line_search(self,f,x,p):
        """
        Perform a line serach of function f from point x along direction p
        """
        a = 1.
        nu = .9
        fx = f.eval(x)
        fx0 = fx

        it = 0
        while it < self.max_lin_step:
            xk = x+a*p
            fxk = f.eval(xk)
            if it > 0:
                if fxk > fx and fx < fx0:
                    #starting to climb again, return the second to last value
                    return a / nu
            fx = fxk
            a *= nu
            if a < self.line_tol:
                break
            it += 1
        return a

