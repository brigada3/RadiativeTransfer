#!/usr/bin/env python3
# -*-coding: utf-8 -*-

import sympy


N = 2001
N1 = 240*N
G = 0.999


def f(s, g=G):
    return g**s


def epsilon_x(s, r):
    return sympy.sqrt(s*(s+2*r))


def kappa_x(s, r, w0):
    math_func = (2*(s+r)+1) * (1-w0*f(s+r))
    return math_func


def fi(v):
    math_func = sympy.Piecewise(
        (sympy.atan(sympy.im(v-1j)/sympy.re(v-1j)), sympy.re(v-1j) > 0),
        (-sympy.pi/2, sympy.And(sympy.Eq(sympy.re(v-1j), 0), sympy.im(v-1j) < 0)),
        (-sympy.pi/2 - sympy.atan(sympy.re(v-1j)/sympy.im(v-1j)), sympy.And(sympy.re(v-1j)<0, sympy.im(v-1j)<0)),
        (-sympy.pi, sympy.And(sympy.re(v-1j)<0, sympy.Eq(sympy.im(v-1j), 0))),
        (-sympy.pi + sympy.atan(sympy.im(v-1j)/sympy.re(v-1j)), sympy.And(sympy.re(v-1j)<0, sympy.im(v-1j)>0))
    )
    return math_func


def fi_1(v):
    math_func = sympy.Piecewise(
        (sympy.atan(sympy.im(v+1j)/sympy.re(v+1j)), sympy.re(v+1j) > 0),
        (sympy.pi/2, sympy.And(sympy.Eq(sympy.re(v+1j), 0), sympy.im(v+1j) > 0)),
        (sympy.pi/2 - sympy.atan(sympy.re(v+1j)/sympy.im(v+1j)), sympy.And(sympy.re(v+1j)<0, sympy.im(v+1j)>0)),
        (sympy.pi, sympy.And(sympy.re(v+1j)<0, sympy.Eq(sympy.im(v+1j), 0))),
        (sympy.pi + sympy.atan(sympy.im(v+1j)/sympy.re(v+1j)), sympy.And(sympy.re(v+1j)<0, sympy.im(v+1j)<0))
    )

    return math_func


def q_x(s, r, w0, n=N):
    math_func = sympy.Piecewise(
        ((s*(s+2*r))/((2*(s+r)+1)*(2*(s+r+1)+1)*(1-w0*f(s+r)*(1-w0*f(s+r+1)))), sympy.And(s+r+1<=n, n>=2)),
        ((s*(s+2*r))/((2*(s+r)+1)*(2*(s+r+1)+1)*(1-w0*f(s+r))), sympy.Eq(s+r, n)),
        ((s*(s+2*r))/((2*(s+r)+1)*(2*(s+r+1)+1)), s >= n-r+1)
    )

    return math_func


def gamma(l, v, r, n=N, n1=N1):
    ind = n1
    result = (
        (-1 +
            (sympy.exp(1)**(fi(v) + fi_1(v))*1j/2) *
            ( (1+sympy.re(v)**2-sympy.im(v)**2)**2 + 4*sympy.re(v)**2*sympy.im(v)**2)**0.5
        ) / 2 + 1
    )

    while ind >= l:
        result = 1 + a(ind, v, r, n)/result
        ind -= 1

    return result


def gamma_x(s, v, r, w0, n=N, n1=N1):
    if s == n-r+1:
        return gamma(1, v, n, r, n, n1)
    elif s.is_Number and r.is_Number:
        gamma_x_s_plus_1 = gamma_x(s+1, v, r, w0, n, n1)
    else:
        gamma_x_s_plus_1 = sympy.Function('gamma_x')(s+1, v, r, w0, n, n1)
    return 1 + q_x(s, r, w0, n)*v/gamma_x_s_plus_1


def h_x(alpha, v, r, w0, n=N, n1=N1):
    math_func = (
        1j*v*epsilon_x(alpha, r) *
        (kappa_x(alpha, r, w0) * gamma_x(alpha, v**2, r, w0, n, n1))**(-1)
    )
    return math_func


def a(l, v, r, n=N):
    math_func = ((n+l)**2-r**2)*v**2 /( (2*(n+l)+1) * (2*(n+l+1)+1) )

    return math_func


def Px(i, r, mu):
    if i == r-1:
        return 0
    elif i == r:
        return (1-mu**2)**(r/2) * (2)**(-r) * sympy.sqrt(sympy.factorial(2*r)/sympy.factorial(r)**2)
    elif i.is_Number:
        Px_i_minus_1 = Px(i-1, r, mu)
        Px_i_minus_2 = Px(i-2, r, mu)
    else:
        Px_i_minus_1 = sympy.Function('Px')(i-1, r, mu)
        Px_i_minus_2 = sympy.Function('Px')(i-2, r, mu)

    return (
        ( 2*i*mu*Px_i_minus_1 - sympy.sqrt((i-r-1)*(i+r-1)) * Px_i_minus_2 ) *
        ( (i-r)*(i+r) )**(-1/2)
    )

def psi_x(j, v, r, w0):
    l = sympy.symbols('l')

    math_func = sympy.Piecewise(
        ( 1, sympy.Eq(j, 0)),
        ( sympy.Product(h_x(l, v, r, w0), (l, 1, j)), True)
    )

    return math_func


def nu_x(alpha, v, r, w0, mu_1, n=N, n1=N1):
    if alpha == n-r+1:
        return 0
    elif alpha == n-r:
        return (
            ((2*n+1)*f(n)*Px(alpha+r, r, mu_1)) /
            (kappa_x(n-r, r, w0)*gamma_x(n-r, v**2, r, w0, n, n1))
        )
    elif alpha.is_Number:    
        nu_x_alpha_plus_1 = nu_x(alpha+1, v, w0, mu_1)
    else:
        nu_x_alpha_plus_1 = sympy.Function('nu_x')(alpha+1, v, w0, mu_1)

    return (
        (kappa_x(alpha, r, w0)*gamma_x(alpha, v**2, r, w0, n, n1))**(-1)*
        (1j*v*epsilon_x(alpha+1, r)*nu_x_alpha_plus_1 + 
            (2*(alpha+r)+1)*f(alpha+r)*Px(alpha+r, r, mu_1)
        )
    )


def b(alpha, v, w0, mu_1, r, n=N, n1=N1):
    if alpha == 0 and n-r >= 0:
        j = sympy.symbols('j')
        return (
            (kappa_x(sympy.Integer('0'), r, w0)*gamma_x(sympy.Integer('0'), v**2, r, w0, n, n1))**(-1)*
            sympy.Sum(
                psi_x(j, v, r, w0)*(2*(j+r)+1)*f(j+r)*Px(j+r, r, mu_1),
                (j, 0, n-r)
            )
        )
    elif alpha.is_Number and 1 <= alpha <= n-r and n-r > 1:
        b_alpha_minus_1 = b(alpha-1, v, w0, mu_1, r, n, n1)
    else:
        b_alpha_minus_1 = sympy.Function('b')(alpha-1, v, w0, mu_1, n, n1)
    
    return (
        h_x(alpha, v, r, w0) * b_alpha_minus_1 +
        nu_x(alpha, v, r, w0, mu_1, n, n1)
    )

def F(s, r, n=N, n1=N1):
    math_func = (0.5)*sympy.Sum(
        (2*(s+r+1))*b(s, v, w0, mu_1, r, n, n1)*Px(s+r, r, mu),
        (s, 0, n)
    )

    return
