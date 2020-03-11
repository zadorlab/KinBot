import sys
import numpy as np


class BFGS:
    """
    Class to execute a BFGS optimization to a local minimum
    """
    def __init__(self, step_tol=1E-7, grad_tol=1E-7, line_tol=1E-10,
                 inhess=None, max_step=100, max_lin_step=1000,
                 use_grad_tol=1, use_step_tol=1):
        """
        Initialize the BFGS algorithm

        step_tol: tolerance on the minimum step size to continue the
        optimization (absolute L2 norm of the step)

        grad_tol: tolerance on the minimum gradient to continue the
        optimization (absolute L2 norm of the gradient)

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
        if inhess is None:
            self.inhess = []
        else:
            self.inhess = inhess
        self.max_step = max_step
        self.max_lin_step = max_lin_step
        self.use_grad_tol = use_grad_tol
        self.use_step_tol = use_step_tol

        if not use_grad_tol and not use_step_tol:
            sys.exit('Cannot execute an optimization if neither the step ' +
                     'nor the gradient tolerance can be used')

    def converged(self, step, grad):
        if np.linalg.norm(step) < self.step_tol and self.use_step_tol:
            return 1
        if np.linalg.norm(grad) < self.grad_tol and self.use_grad_tol:
            return 1
        return 0

    def optimize(self, f, x):
        """
        Optimize a function to the closest local minimum using the BFGS method.

        f: a function of x which needs to have
            f.eval(x): function evaluation of f
            f.gradient(x): gradient df/dx in x

        x: initial guess
        """
        # intermediate coordinates
        x_i = [x]

        if len(self.inhess) == 0:
            H = np.eye(len(x))
        else:
            H = self.inhess

        Hinv = np.linalg.inv(H)

        g = f.gradient(x)
        # intermediate forces
        g_i = [g]
        p = np.dot(np.linalg.inv(H), -g)

        it = 0
        while it < self.max_step:
            ak = self.line_search(f, x, p)
            sk = ak * p
            xk = x + sk
            gk = f.gradient(xk)
            yk = gk - g

            sy = np.dot(sk, yk)
            if sy:
                t1 = (sy + np.dot(yk, np.dot(Hinv, yk))) * np.outer(sk, sk)
                t1 = t1 / np.power(sy, 2)
                t2 = np.dot(Hinv, np.outer(yk, sk))
                t2 = t2 + np.dot(np.outer(sk, yk), Hinv)
                t2 = t2 / sy
                Hinvk = Hinv + t1 - t2
            else:
                Hinvk = Hinv

            x_i.append(xk)
            g_i.append(gk)

            converged = self.converged(sk, gk)
            if converged:
                return xk, x_i, g_i
            g = gk
            x = xk
            p = np.dot(Hinvk, -gk)
            Hinv = Hinvk

            it += 1
        # return final x values, list of geometric steps
        # and list of gradients
        return xk, x_i, g_i

    def line_search(self, f, x, p):
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
                    # starting to climb again, return the second to last value
                    return a / nu
            fx = fxk
            a *= nu
            if a < self.line_tol:
                break
            it += 1
        return a
