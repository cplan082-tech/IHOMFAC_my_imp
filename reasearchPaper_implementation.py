# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 18:07:01 2022

@author: clive
"""
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
import class_IHOMFAC
import pandas as pd
# plt.rcParams['figure.figsize'] = [20, 10]

def desired_output(t):
    if (t <= 300) or (t > 700):
        y = 0.5*(-1)**(round(t/100))

    else:
        y = 0.5*np.sin(t*pi/100) + 0.3*np.cos(t*pi/50)

    return y

sim_len = 1000
t = [i for i in range(sim_len + 1)]
lst_y_setpoint = [desired_output(t_i) for t_i in t]
# lst_y_setpoint = [0.5*(-1)**(round(t[i]/200)) for i in range(len(t))]

lst_y = [0]*3 # list of measured outputs
lst_u = [0]*2 # list of measured control inputs

eta = 0.8
lam = 0.1
mu = 0.01
ro = 0.8
eps = 1e-15
alpha = [0.5, 0.25, 1/8, 1/8]
beta = [0.5, 0.25, 1/8, 1/16, 1/32, 1/32]

obj_IHOMFAC = class_IHOMFAC.IHOMFAC(eta, lam, mu, ro, eps, alpha, beta, phi_init=0.1) # best phi_init = 0.1

y_meas = 0 # initial output measurment

for k in t:

    obj_IHOMFAC.run_algo(y_meas, lst_y_setpoint[k])
    
    lst_u.append(obj_IHOMFAC.u[0].reshape(())) # TODO: check if this is giving what we want
     
      # lst index [-1] signifies time t, [-2] signifies time (t-1) and so on
    if k <= 500:
        y_meas = lst_y[-1]/(1 + lst_y[-1]**2) + lst_u[-1]**3
       
    else:
        y_meas = ((lst_y[-1]*lst_y[-2]*lst_y[-3]*lst_u[-2]*
                (lst_y[-3] - 1) + round(k/500)*lst_u[-1])/
                (1 + lst_y[-2]**2 + lst_y[-3]**2))
   
    # append new output measurment to end of list
    lst_y.append(y_meas)


[lst_u.pop(0) for _ in range(2)] # removes first two elements from lst_u (init vals)
[lst_y.pop(0) for _ in range(3)] # removes first three elements from lst_y (init vals)


plt.plot(t, lst_y_setpoint)
plt.plot(t, lst_y)
plt.grid()
plt.ylim((-2,2))
plt.xlim((0,sim_len))
plt.show()

# plt.plot(t, lst_u[2:])
# plt.show()

# zoom_range_low = 130
# zoom_range_high = 150 # 252
# plt.plot(t[zoom_range_low:zoom_range_high], lst_y_setpoint[zoom_range_low:zoom_range_high])
# plt.plot(t[zoom_range_low:zoom_range_high], lst_y[(zoom_range_low ):(zoom_range_high)])


