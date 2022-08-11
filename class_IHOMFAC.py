# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 18:07:27 2022

@author: clive
"""

import numpy as np
from scipy.ndimage import shift

class IHOMFAC:
    def __init__(self, eta, lam, mu, ro, eps, alpha, beta, phi_init=0, y_init=0, u_init=0):
        self.eta = eta
        self.lam = lam
        self.mu = mu
        self.ro = ro
        self.eps = eps
        
        self.alpha = np.array(alpha).reshape(1,-1) # row vec
        self.beta = np.array(beta).reshape(1,-1) # row vec
        
        self.y = np.array([y_init]*2)
        self.u = np.array([u_init]*len(alpha)).reshape(-1,1) # col vec
        self.u[1] = 1 # for testing
        
        self.phi_init = phi_init
        self.phi_init_sign = np.sign(phi_init)
        self.phi = np.array([phi_init]*len(beta)).reshape(-1,1) # col vec
        
        
    def param_est(self, y_meas):
        self.y = shift(self.y, 1, cval=y_meas) # shifts y arr right and replaces first val with new y meas
        
        # TESTINGVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
#         print(np.shape(self.beta))
#         print(np.shape(self.phi))
        sum_beta = self.beta@self.phi # TODO: Reshape this to ()
        del_u_t_1 = self.u[1] - self.u[2]
        del_y_t = self.y[0] - self.y[1] # TODO: Reshape this to ()
        
        # New phi estimate
        phi_est = (sum_beta + (self.eta*del_u_t_1)*(del_y_t - del_u_t_1*sum_beta)/(self.mu + del_u_t_1**2))
        
        if (phi_est <= self.eps) or (abs(del_u_t_1) <= self.eps) or (np.sign(phi_est) != self.phi_init_sign):
            phi_est = self.phi_init
            
        self.phi = shift(self.phi, [1,0], cval=phi_est) # shifts phi arr right and replaces first val with new phi est
        
        
    def control_input(self, y_setpoint):
        sum_alpha = self.alpha@self.u
        coeff_denum = self.lam + self.phi[0]**2
        error = y_setpoint - self.y[0]
        
        u_calc = ((self.phi[0]**2)*self.u[1] + self.lam*sum_alpha + self.ro*self.phi[0]*error)/coeff_denum
        self.u = shift(self.u, [1,0], cval=u_calc) # shifts u arr right and replaces first val with new control input u
        
        
    def run_algo(self, y_meas, y_setpoint):
        self.param_est(y_meas)
        self.control_input(y_setpoint)
        
        
        