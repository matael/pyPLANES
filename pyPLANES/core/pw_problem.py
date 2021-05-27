#! /usr/bin/env python
# -*- coding:utf8 -*-
#
# pw_classes.py
#
# This file is part of pyplanes, a software distributed under the MIT license.
# For any question, please contact one of the authors cited below.
#
# Copyright (c) 2020
# 	Olivier Dazel <olivier.dazel@univ-lemans.fr>
# 	Mathieu Gaborit <gaborit@kth.se>
# 	Peter Göransson <pege@kth.se>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
from termcolor import colored
import numpy as np
import numpy.linalg as LA

from numpy import pi

import matplotlib.pyplot as plt

from mediapack import Air, Fluid

# from pyPLANES.utils.io import initialisation_out_files_plain
from pyPLANES.core.calculus import Calculus
from pyPLANES.pw.multilayer import MultiLayer

from pyPLANES.pw.pw_layers import *
from pyPLANES.pw.pw_interfaces import *

class PwProblem(Calculus, MultiLayer):
    """
        Plane Wave Problem 
    """ 
    def __init__(self, **kwargs):
        Calculus.__init__(self, **kwargs)
        self.theta_d = kwargs.get("theta_d", 0.0)
        self.method = kwargs.get("method", False)
        if self.method.lower() in ["recursive", "jap", "recursive method"]:
            self.method = "Recursive Method"
            self.info_file.write("Plane Wave solver // Recursive method\n")
            if self.theta_d == 0:
                self.theta_d = 1e-12 
        else: 
            self.method = "Global Method"
            self.info_file.write("Plane Wave solver // Global method\n")
        # Out files
        if self.method == "Global Method":
            self.out_file_extension = "GM"
        elif self.method == "Recursive Method":
            self.out_file_extension = "RM"

        assert "ml" in kwargs
        ml = kwargs.get("ml")

        MultiLayer.__init__(self, ml)
        self.termination = kwargs.get("termination", "rigid")

        self.add_excitation_and_termination(self.method, self.termination)

        # Calculus variable (for pylint)
        self.kx, self.ky, self.k = None, None, None
        self.R, self.T = None, None


    def update_frequency(self, omega):
        Calculus.update_frequency(self, omega)
        self.kx = np.array([omega*np.sin(self.theta_d*np.pi/180)/Air.c])
        self.ky = np.array([omega*np.cos(self.theta_d*np.pi/180)/Air.c])
        self.k_air = omega/Air.c
        MultiLayer.update_frequency(self, omega, self.kx)

    def create_linear_system(self, omega):
        Calculus.create_linear_system(self, omega)
        if self.method == "Recursive Method":
            if self.termination == "transmission":
                self.Omega, self.back_prop = self.interfaces[-1].Omega()
                # print(self.interfaces[-1])
                # print(self.Omega)
                for i, _l in enumerate(self.layers[::-1]):
                    next_interface = self.interfaces[-i-2]
                    # print(_l)
                    # print(next_interface)
                    _l.Omega_plus, _l.Xi = _l.transfert(self.Omega)
                    self.back_prop = self.back_prop@_l.Xi
                    self.Omega, next_interface.Tau = next_interface.transfert(_l.Omega_plus)
                    self.back_prop = self.back_prop@next_interface.Tau
            else: # Rigid backing
                self.Omega = self.interfaces[-1].Omega()
                for i, _l in enumerate(self.layers[::-1]):
                    next_interface = self.interfaces[-i-2]
                    _l.Omega_plus, _l.Xi = _l.transfert(self.Omega)
                    self.Omega, next_interface.Tau = next_interface.transfert(_l.Omega_plus)

        elif self.method == "Global Method":
            self.A = np.zeros((self.nb_PW-1, self.nb_PW),dtype=complex)
            i_eq = 0
            # Loop on the interfaces
            for _int in self.interfaces:
                if self.method == "Global Method":
                    i_eq = _int.update_M_global(self.A,i_eq)
            self.F = -self.A[:, 0]*np.exp(1j*self.ky*self.layers[0].d) # - is for transposition, exponential term is for the phase shift
            self.A = np.delete(self.A, 0, axis=1)

    def solve(self):
        Calculus.solve(self)
        if self.method == "Recursive Method":
            self.Omega = self.Omega.reshape(2)
            _ = 1j*(self.ky[0]/self.k_air)/(2*pi*self.f*Air.Z)
            det = -self.Omega[0]+_*self.Omega[1]
            self.R = (self.Omega[0]+_*self.Omega[1])/det
            self.X_0_minus = 2*_/det
            if self.termination == "transmission":
                Omega_end = (self.back_prop*self.X_0_minus).flatten()
                self.T = Omega_end[0]
        elif self.method == "Global Method":
            self.X = LA.solve(self.A, self.F)
            self.R = self.X[0]
            if self.termination == "transmission":
                self.T = self.X[-1]
            else:
                self.T = None
        if self.print_result:
            _text = "R={:+.15f}".format(self.R)
            if self.termination == "transmission":
                _text += " / T={:+.15f}".format(self.T)
            print(_text)

    def plot_solution(self):
        if self.method == "Global Method":
            for _l in self.layers[1:]:
                _l.plot_solution_global(self.plot,self.X[_l.dofs-1])  
        else:
            x = np.array([self.X_0_minus]) # Information vector at incident interface  x^-
            for i, _l in enumerate(self.layers):
                x = self.interfaces[i].Tau @ x # Transfert through the interface x^+
                q = LA.solve(_l.SV, _l.Omega_plus@x)
                _l.plot_solution_recursive(self.plot, q)
                x = _l.Xi@x # Transfert through the layer x^-_{+1}

    def write_out_files(self):
        self.out_file.write("{:.12e}\t".format(self.f))
        self.out_file.write("{:.12e}\t".format(self.R.real))
        self.out_file.write("{:.12e}\t".format(self.R.imag))
        if self.termination == "transmission":
            self.out_file.write("{:.12e}\t".format(self.T.real))
            self.out_file.write("{:.12e}\t".format(self.T.imag))
        self.out_file.write("\n")

    def load_results(self):

        name_file_with_method = self.out_file_name.split(".")
        name_file_with_method.insert(1, self.out_file_extension)
        name_file_with_method = ".".join(name_file_with_method)

        data = np.loadtxt(name_file_with_method)
        f  = data[:, 0]
        R  = data[:, 1] + 1j*data[:, 2]
        if self.termination == "transmission":
            T = data[:,3] + 1j*data[:,4]
        else:
            T= False
        return f, R, T

    def compute_error(self, f_ref, R_ref, T_ref, plot_RT=False, eps=1e-8):
        # self.frequencies = f_ref
        self.resolution()
        R, T = self.load_results()[1:] 
        error_R = LA.norm(R-R_ref)
        error_T = LA.norm(T-T_ref)
        _ = np.sqrt(error_R**2+error_T**2)/len(f_ref)
        if plot_RT:
            self.plot_results()

        if _< eps:
            print("Overall error = {}".format(_) + "\t"*6 + "["+ colored("OK", "green")  +"]")
        else:
            print("Overall error = {}".format(_) + "\t"*6 + "["+ colored("Fail", "red")  +"]")
        return error_R, error_T

    def plot_results(self):
        if self.method == "Recursive Method":
            method = " RM"
            marker = "."
        elif self.method == "Global Method":
            method = " RM"
            marker = "+"
        f, R, T = self.load_results()
        plt.figure(self.name_project + "/ Reflection coefficient")
        plt.plot(f, np.real(R),"r"+marker,label="Re(R)"+method)
        plt.plot(f, np.imag(R),"b"+marker,label="Im(R)"+method)
        plt.legend()
        if self.termination == "transmission":
            plt.figure(self.name_project + "/ Transmission coefficient")
            plt.plot(f, np.real(T),"r"+marker,label="Re(T)"+method)
            plt.plot(f, np.imag(T),"b"+marker,label="Im(T)"+method)
            plt.legend()