#! /usr/bin/env python
# -*- coding:utf8 -*-
#
# problem.py
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


import numpy as np
from mediapack import Air


class Calculus():
    """
    Calculus - Base class for all computations

    This class is mostly storing important information about ongoing
    computations and is intended to be used as a base for other computation
    classes

    Parameters
    ----------
    name_project : str
        Name of the project
    frequencies : np.ndarray
        Vectors of frequencies to compute
    theta_d : float, None
        Angle of incidence (0 if ignored)
    outfiles_directory : str, False
        Specific directory in which to create the output files
    plot_results : list, False
        Controls which plots to create
        todo: define the mapping
    """

    def __init__(self, name_project="unnamed_project", frequencies=np.array([440]),
                 theta_d=False, outfiles_directory=False, plot_results=False, materials=None):
        self.frequencies = self.init_vec_frequencies(frequencies)
        self.current_frequency = None
        self.omega = None
        self.theta_d = theta_d
        self.name_project = name_project
        self.outfiles_directory = outfiles_directory
        self.plot = plot_results
        self.materials = materials if materials is not None else {}

    def init_vec_frequencies(self, frequency):
        if frequency[2] > 0:
                frequencies = np.linspace(frequency[0], frequency[1], frequency[2])
        elif frequency[2]<0:
            frequencies = np.logspace(np.log10(frequency[0]),np.log10(frequency[1]),abs(frequency[2]))
        # else % Case of complex frequency
        #     temp_1=linspace(frequency.min,frequency.max,frequency.nb(1));
        #     temp_2=linspace(frequency.min_imag,frequency.max_imag,frequency.nb(2));
        #     frequency.vec=[];
        #     for ii=1:frequency.nb(2)
        #         frequency.vec=[frequency.vec temp_1+1j*temp_2(ii)];
        #     end
        #     frequency.nb=frequency.nb(1)*frequency.nb(2);
        return frequencies

    def update_frequency(self, f):
        self.current_frequency = f
        self.omega = 2*np.pi*f


class FemCalculus(Calculus):

    def __init__(self, **kwargs):
        Calculus.__init__(self, **kwargs)
        self.out_file = self.name_project + ".FEM.txt"
        self.info_file = self.name_project + ".info.FEM.txt"
        self.F_i, self.F_v = None, None
        self.A_i, self.A_j, self.A_v = None, None, None
        self.A_i_c, self.A_j_c, self.A_v_c = None, None, None
        self.T_i, self.T_j, self.T_v = None, None, None
        self.modulus_reflex, self.modulus_trans, self.abs = None, None, None

    def update_frequency(self, f):
        Calculus.update_frequency(self, f)
        self.F_i, self.F_v = [], []
        self.A_i, self.A_j, self.A_v = [], [], []
        self.A_i_c, self.A_j_c, self.A_v_c = [], [], []
        self.T_i, self.T_j, self.T_v = [], [], []
        self.kx = (self.omega/Air.c)*np.sin(self.theta_d*np.pi/180)
        self.ky = (self.omega/Air.c)*np.cos(self.theta_d*np.pi/180)
        self.delta_periodicity = np.exp(-1j*self.kx*self.period)
        self.nb_dofs = self.nb_dof_FEM
        for _ent in self.model_entities:
            _ent.update_frequency(self.omega)
        self.modulus_reflex, self.modulus_trans, self.abs = 0, 0, 1


class PwCalculus(Calculus):
    def __init__(self, **kwargs):
        Calculus.__init__(self, **kwargs)
        self.out_file = self.name_project + ".PW.txt"
        self.info_file = self.name_project + ".info.PW.txt"

    def update_frequency(self, f):
        Calculus.update_frequency(self, f)
        for _l in self.layers:
            _l.medium.update_frequency(self.omega)
        self.kx = self.omega*np.sin(self.theta_d*np.pi/180)/Air.c
        self.ky = self.omega*np.cos(self.theta_d*np.pi/180)/Air.c
        self.k = self.omega/Air.c
