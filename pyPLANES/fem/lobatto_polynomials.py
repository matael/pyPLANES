#! /usr/bin/env python
# -*- coding:utf8 -*-
#
# lobatto_polynomials.py
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

def lobatto(order,x):
    """Lobatto shape functions

    Eq (1.51) of Solin

    Parameters
    ----------
    order : integer
        degree of the polynomial
    x : ndarray
        position of the Gauss points

    Returns
    -------
    p_lobatto : ndarray
        value of the Lobatto polynomial at x 
    d_lobatto : ndarray
        value of the Lobatto polynomial first derivative at x      
    """
    if order == 0:
        p_lobatto = 1-x
        t_lobatto = np.array([1, -1], dtype=np.float32)
        d_lobatto = -x**0
        p_lobatto /= 2
        d_lobatto /= 2
        t_lobatto /= 2
    elif order == 1:
        p_lobatto = 1+x
        t_lobatto = np.array([1, 1], dtype=np.float32)
        d_lobatto = x**0
        p_lobatto /= 2
        d_lobatto /= 2
        t_lobatto /= 2
    elif order == 2:
        p_lobatto = -1+x**2
        t_lobatto = np.array([-1, 0, 1], dtype=np.float32)
        d_lobatto = 2*x
        p_lobatto *= np.sqrt(3./2.)/2.
        d_lobatto *= np.sqrt(3./2.)/2.
        t_lobatto *= np.sqrt(3./2.)/2.
    elif order == 3:
        p_lobatto = -x+x**3
        t_lobatto = np.array([0, -1, 0, 1], dtype=np.float32)
        d_lobatto = -1+3*x**2
        p_lobatto *= np.sqrt(5./2.)/2.
        t_lobatto *= np.sqrt(5./2.)/2.
        d_lobatto *= np.sqrt(5./2.)/2.
    elif order == 4:
        p_lobatto = 1+x**2*(-6+5*x**2)
        t_lobatto = np.array([1, 0, -6, 0, 5], dtype=np.float32)
        d_lobatto = x*(-12+20*x**2)
        p_lobatto *= np.sqrt(7./2.)/8.
        t_lobatto *= np.sqrt(7./2.)/8.
        d_lobatto *= np.sqrt(7./2.)/8.
    elif order == 5:
        p_lobatto = x*(3+x**2*(-10+7*x**2))
        t_lobatto = np.array([0, 3, 0, -10, 0, 7], dtype=np.float32)
        d_lobatto = 3+x**2*(-30+35*x**2)
        p_lobatto *= (3/8)/np.sqrt(2)
        t_lobatto *= (3/8)/np.sqrt(2)
        d_lobatto *= (3/8)/np.sqrt(2)
    elif order == 6:
        p_lobatto = (-1)+x**2*(15+x**2*((-35)+21*x**2))
        t_lobatto = np.array([-1, 0, 15, 0, -35, 0, 21], dtype=np.float32)
        d_lobatto = x*(30+x**2*((-140)+126*x**2))
        p_lobatto *= (1/16)*np.sqrt(11/2)
        t_lobatto *= (1/16)*np.sqrt(11/2)
        d_lobatto *= (1/16)*np.sqrt(11/2)
    elif order == 7:
        p_lobatto = x*((-5)+x**2*(35+x**2*((-63)+33*x**2)))
        t_lobatto = np.array([0, -5, 0, 35, 0, -63, 0, 33], dtype=np.float32)
        d_lobatto = (-5)+x**2*(105+x**2*((-315)+231*x**2))
        p_lobatto *= (1/16)*np.sqrt(13/2)
        t_lobatto *= (1/16)*np.sqrt(13/2)
        d_lobatto *= (1/16)*np.sqrt(13/2)
    elif order == 8:
        p_lobatto = 5+x**2*((-140)+x**2*(630+x**2*((-924)+429*x**2)))
        t_lobatto = np.array([5, 0, -140, 0, 630, 0, -924, 0, 429], dtype=np.float32)
        d_lobatto = x*((-280)+x**2*(2520+x**2*((-5544)+3432*x**2)))
        p_lobatto *= (1/128)*np.sqrt(15/2)
        t_lobatto *= (1/128)*np.sqrt(15/2)
        d_lobatto *= (1/128)*np.sqrt(15/2)
    elif order == 9:
        p_lobatto = x*(35+x**2*((-420)+x**2*(1386+x**2*((-1716)+715*x**2))))
        t_lobatto = np.array([0, 35, 0, -420, 0, 1386, 0, -1716, 0, 715], dtype=np.float32)
        d_lobatto = 35+x**2*((-1260)+x**2*(6930+x**2*((-12012)+6435*x**2)))
        p_lobatto *= (1/128)*np.sqrt(17/2)
        t_lobatto *= (1/128)*np.sqrt(17/2)
        d_lobatto *= (1/128)*np.sqrt(17/2)
    elif order == 10:
        p_lobatto = (-7)+x**2*(315+x**2*((-2310)+x**2*(6006+x**2*((-6435)+2431*x**2))))
        t_lobatto = np.array([-7, 0, 315, 0, -2310, 0, 6006, 0, -6435, 0, 2431], dtype=np.float32)
        d_lobatto = x*(630+x**2*((-9240)+x**2*(36036+x**2*((-51480)+24310*x**2))))
        p_lobatto *= (1/256)*np.sqrt(19/2)
        t_lobatto *= (1/256)*np.sqrt(19/2)
        d_lobatto *= (1/256)*np.sqrt(19/2)
    elif order == 11:
        p_lobatto = x*((-63)+x**2*(1155+x**2*((-6006)+x**2*(12870+x**2*((-12155)+4199*x**2)))))
        t_lobatto = np.array([0, -63, 0, 1155, 0, -6006, 0, 12870, 0, -12155, 0, 4199], dtype=np.float32)
        d_lobatto = (-63)+x**2*(3465+x**2*((-30030)+x**2*(90090+x**2*((-109395)+46189*x**2))))
        p_lobatto *= (1/256)*np.sqrt(21/2)
        t_lobatto *= (1/256)*np.sqrt(21/2)
        d_lobatto *= (1/256)*np.sqrt(21/2)
    elif order == 12:
        p_lobatto = 21+x**2*((-1386)+x**2*(15015+x**2*((-60060)+x**2*(109395+x**2*((-92378)+29393*x**2)))))
        t_lobatto = np.array([21, 0, -1386, 0, 15015, 0, -60060, 0, 109395, 0, -92378, 0, 29393], dtype=np.float32)
        d_lobatto = x*((-2772)+x**2*(60060+x**2*((-360360)+x**2*(875160+x**2*((-923780)+352716*x**2)))))
        p_lobatto *= (1/1024)*np.sqrt(23/2)
        t_lobatto *= (1/1024)*np.sqrt(23/2)
        d_lobatto *= (1/1024)*np.sqrt(23/2)
    elif order == 13:
        p_lobatto = x*(231+x**2*((-6006)+x**2*(45045+x**2*((-145860)+x**2*(230945+x**2*((-176358)+52003*x**2))))))
        t_lobatto = np.array([0, 231, 0, -6006, 0, 45045, 0, -145860, 0, 230945, 0, -176358, 0, 52003], dtype=np.float32)
        d_lobatto = 231+x**2*((-18018)+x**2*(225225+x**2*((-1021020)+x**2*(2078505+x**2*((-1939938)+676039*x**2)))))
        p_lobatto *= (5/1024)*np.sqrt(2)**(-1/2)
        t_lobatto *= (5/1024)*np.sqrt(2)**(-1/2)
        d_lobatto *= (5/1024)*np.sqrt(2)**(-1/2)
    elif order == 14:
        p_lobatto = (-33)+x**2*(3003+x**2*((-45045)+x**2*(255255+x**2*((-692835)+x**2*(969969+x**2*((-676039)+185725*x**2))))))
        t_lobatto = np.array([-33, 0, 3003, 0, -45045, 0, 255255, 0, -692835, 0, 969969, 0, -676039, 0, 185725], dtype=np.float32)
        d_lobatto = x*(6006+x**2*((-180180)+x**2*(1531530+x**2*((-5542680)+x**2*(9699690+x**2*((-8112468)+2600150*x**2))))))
        p_lobatto *= (3/2048)*np.sqrt(3/2)
        t_lobatto *= (3/2048)*np.sqrt(3/2)
        d_lobatto *= (3/2048)*np.sqrt(3/2)
    elif order == 15:
        p_lobatto = x*((-429)+x**2*(15015+x**2*((-153153)+x**2*(692835+x**2*((-1616615)+x**2*(2028117+x**2*((-1300075)+334305*x**2)))))))
        t_lobatto = np.array([0, -429, 0, 15015, 0, -153153, 0, 692835, 0, -1616615, 0, 2028117, 0, -1300075, 0, 334305], dtype=np.float32)
        d_lobatto = (-429)+x**2*(45045+x**2*((-765765)+x**2*(4849845+x**2*((-14549535)+x**2*(22309287+x**2*((-16900975)+5014575*x**2))))))
        p_lobatto *= (1/2048)*np.sqrt(29/2)
        t_lobatto *= (1/2048)*np.sqrt(29/2)
        d_lobatto *= (1/2048)*np.sqrt(29/2)
    return p_lobatto, d_lobatto, t_lobatto

def lobatto_kernels(order, x):
    """Lobatto kernels

    Equation (1.53) of Solin

    Parameters
    ----------
    order : integer
        degree of the polynomial
    x : ndarray
        position of the Gauss points
    Returns
    -------
    phi : ndarray
        value of the Lobatto kernel at x 
    dphi : ndarray
        value of the Lobatto kernel first derivative at x      
    """
    if order == 0:
        phi = -np.sqrt(6)*x**0
        dphi = 0*x
    elif order == 1:
        phi = -np.sqrt(10.)*x
        dphi = -np.sqrt(10.)*x**0
    elif order == 2:
        phi = -1+5*x**2
        dphi = 10*x
        phi *= -np.sqrt(7./2.)/2.
        dphi *= -np.sqrt(7./2.)/2.
    elif order == 3:
        phi = x*(3-7*x**2)
        dphi = 3-21*x**2
        phi *= (3/2)/np.sqrt(2)
        dphi *= (3/2)/np.sqrt(2)
    elif order == 4:
        phi = -1+x**2*(14-21*x**2)
        dphi = x*(28+(-84)*x**2)
        phi *= (1/4)*np.sqrt(11./2.)
        dphi *= (1/4)*np.sqrt(11./2.)
    elif order == 5:
        phi = x*(-5+x**2*(30+-33*x**2))
        dphi = (-5)+x**2*(90+(-165)*x**2)
        phi *= np.sqrt(13/2)/4
        dphi *= np.sqrt(13/2)/4
    elif order == 6:
        phi = 5+x**2*((-135)+x**2*(495+(-429)*x**2))
        dphi = x*((-270)+x**2*(1980+(-2574)*x**2))
        phi = (1/32)*(15/2)**(1/2)*phi
        dphi = (1/32)*(15/2)**(1/2)*dphi
    elif order == 7:
        phi = x*(35+x**2*((-385)+x**2*(1001+(-715)*x**2)))
        dphi = 35+x**2*((-1155)+x**2*(5005+(-5005)*x**2))
        phi = (1/32)*(17/2)**(1/2)*phi
        dphi = (1/32)*(17/2)**(1/2)*dphi
    elif order == 8:
        phi = (-7)+x**2*(308+x**2*((-2002)+x**2*(4004+(-2431)*x**2)))
        dphi = x*(616+x**2*((-8008)+x**2*(24024+(-19448)*x**2)))
        phi = (1/64)*(19/2)**(1/2)*phi
        dphi = (1/64)*(19/2)**(1/2)*dphi
    elif order == 9:
        phi = x*((-63)+x**2*(1092+x**2*((-4914)+x**2*(7956+(-4199)*x**2))))
        dphi = (-63)+x**2*(3276+x**2*((-24570)+x**2*(55692+(-37791)*x**2)))
        phi = (1/64)*(21/2)**(1/2)*phi
        dphi = (1/64)*(21/2)**(1/2)*dphi
    elif order == 10:
        phi = 21+x**2*((-1365)+x**2*(13650+x**2*((-46410)+x**2*(62985+(-29393)*x**2))))
        dphi = x*((-2730)+x**2*(54600+x**2*((-278460)+x**2*(503880+(-293930)*x**2))))
        phi = (1/256)*(23/2)**(1/2)*phi
        dphi = (1/256)*(23/2)**(1/2)*dphi
    elif order == 11:
        phi = x*(231+x**2*((-5775)+x**2*(39270+x**2*((-106590)+x**2*( 124355+(-52003)*x**2)))))
        dphi = 231+x**2*((-17325)+x**2*(196350+x**2*((-746130)+x**2*(1119195+(-572033)*x**2))))
        phi = (5/256)*2**(-1/2)*phi
        dphi = (5/256)*2**(-1/2)*dphi
    elif order == 12:
        phi = (-33)+x**2*(2970+x**2*((-42075)+x**2*(213180+x**2*((-479655)+x**2*(490314+(-185725)*x**2)))))
        dphi = x*(5940+x**2*((-168300)+x**2*(1279080+x**2*((-3837240)+x**2*(4903140+(-2228700)*x**2)))))
        phi = (3/512)*(3/2)**(1/2)*phi
        dphi = (3/512)*(3/2)**(1/2)*dphi
    elif order == 13:
        phi = x*((-429)+x**2*(14586+x**2*((-138567)+x**2*(554268+x**2*((-1062347)+x**2*(965770+(-334305)*x**2))))))
        dphi = (-429)+x**2*(43758+x**2*((-692835)+x**2*(3879876+x**2*((-9561123)+x**2*(10623470+(-4345965)*x**2)))))
        phi = (1/512)*(29/2)**(1/2)*phi
        dphi = (1/512)*(29/2)**(1/2)*dphi
    return phi, dphi

def lobatto_legendre(order, x):
    from scipy.special import eval_legendre
    return order*(eval_legendre(order-1,x)-x*eval_legendre(order,x))

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    xi = np.linspace(-1,1,500)
    for order in range(2,6):
        plt.figure(order)
        plt.plot(xi, lobatto(order,xi)[0])
        plt.plot(xi, lobatto(0,xi)[0]*lobatto(1,xi)[0]*lobatto_kernels(order-2, xi)[0],"r.")
        # plt.plot(xi, lobatto_legendre(order-1, xi),"m+")
        # plt.figure(order*100)
        # plt.plot(xi, lobatto(order,xi)[1])
        # __  = lobatto(0,xi)[1]*lobatto(1,xi)[0]*lobatto_kernels(order-2, xi)[0]
        # __ += lobatto(0,xi)[0]*lobatto(1,xi)[1]*lobatto_kernels(order-2, xi)[0]
        # __ += lobatto(0,xi)[0]*lobatto(1,xi)[0]*lobatto_kernels(order-2, xi)[1]
        # plt.plot(xi, __,"r.")
    plt.show()

