# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 18:07:01 2022

@author: clive
"""
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
import class_IHOMFAC
# from class_IHOMFAC import IHOMFAC
# plt.rcParams['figure.figsize'] = [20, 10]


class data_gen:
    def __init__(self):
        # end of lists below (y and u) represent newest elements while first element is oldest
        self.lst_y = [0]*3 # list of measured outputs
        self.lst_u = [0]*2 # list of measured control inputs
        self.t = 0
        
    def output_meas(self, u):
        # append new control input to end of list
        self.lst_u.append(u)
        
        # lst index [-1] signifies time t, [-2] signifies time (t-1) and so on
        if self.t <= 500:
            y_meas = self.lst_y[-1]/(1 + self.lst_y[-1]**2) + self.lst_u[-1]**3
            
        else:
            y_meas = ((self.lst_y[-1]*self.lst_y[-2]*self.lst_y[-3]*self.lst_u[-2]*
                     (self.lst_y[-3] - 1) + round(self.t/500)*self.lst_u[-1])/
                     (1 + self.lst_y[-2]**2 + self.lst_y[-3]**2))
        
        # append new output measurment to end of list
        self.lst_y.append(y_meas)
        self.t += 1
        return y_meas
    
    

def desired_output(t):
    if (t <= 300) or (t > 700):
        y = 0.5*(-1)**(round(t/100))

    else:
        y = 0.5*np.sin(t*pi/100) + 0.3*np.cos(t*pi/50)

    return y


t = [i for i in range(1001)]
lst_y_setpoint = [desired_output(t_i) for t_i in t]

eta = 0.8
lam = 0.1
mu = 0.01
ro = 0.8
eps = 5e-16
alpha = [0.5, 0.25, 1/8, 1/8]
beta = [0.5, 0.25, 1/8, 1/16, 1/32, 1/32]

obj_data = data_gen()
obj_IHOMFAC = class_IHOMFAC.IHOMFAC(eta, lam, mu, ro, eps, alpha, beta, phi_init=1)

for k in t:
    y_meas = obj_data.output_meas(obj_IHOMFAC.u[0].reshape(()))
    obj_IHOMFAC.run_algo(y_meas, lst_y_setpoint[k])

plt.plot(t, lst_y_setpoint)
plt.plot(t, obj_data.lst_y[3:])
plt.grid()
plt.ylim((-2,2))
plt.xlim((0,1000))