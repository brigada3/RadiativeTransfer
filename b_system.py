#!/usr/bin/env python3
# -*-coding: utf-8 -*-

import sympy


N = 2001
N1 = 240*N
G = 0.999


s, r, w0, i, n = sympy.symbols('s r w0 i n')
alpha, l, j, v = sympy.symbols('alpha l j v')
mu, mu_1 = sympy.symbols('mu mu_1')


f = G**i
epsilon_x = sympy.sqrt(s*(s+2*r))
kappa_x = (2*(s+r)+1)*(1-w0*f.subs(i, s+r))

fi = sympy.Piecewise(
    (sympy.atan(sympy.im(v-1j)/sympy.re(v-1j)), sympy.re(v-1j) > 0),
    (-sympy.pi/2, sympy.And(sympy.Eq(sympy.re(v-1j), 0), sympy.im(v-1j) < 0)),
    (-sympy.pi/2 - sympy.atan(sympy.re(v-1j)/sympy.im(v-1j)), sympy.And(sympy.re(v-1j)<0, sympy.im(v-1j)<0)),
    (-sympy.pi, sympy.And(sympy.re(v-1j)<0, sympy.Eq(sympy.im(v-1j), 0))),
    (-sympy.pi + sympy.atan(sympy.im(v-1j)/sympy.re(v-1j)), sympy.And(sympy.re(v-1j)<0, sympy.im(v-1j)>0))
)

fi_1 = sympy.Piecewise(
    (sympy.atan(sympy.im(v+1j)/sympy.re(v+1j)), sympy.re(v+1j) > 0),
    (sympy.pi/2, sympy.And(sympy.Eq(sympy.re(v+1j), 0), sympy.im(v+1j) > 0)),
    (sympy.pi/2 - sympy.atan(sympy.re(v+1j)/sympy.im(v+1j)), sympy.And(sympy.re(v+1j)<0, sympy.im(v+1j)>0)),
    (sympy.pi, sympy.And(sympy.re(v+1j)<0, sympy.Eq(sympy.im(v+1j), 0))),
    (sympy.pi + sympy.atan(sympy.im(v+1j)/sympy.re(v+1j)), sympy.And(sympy.re(v+1j)<0, sympy.im(v+1j)<0))
)

q_x = sympy.Piecewise(
    ((s*(s+2*r))/((2*(s+r)+1)*(2*(s+r+1)+1)*(1-w0*f.subs(i, s+r)*(1-w0*f.subs(i, s+r+1)))), sympy.And(s+r+1<=N, N>=2)),
    ((s*(s+2*r))/((2*(s+r)+1)*(2*(s+r+1)+1)*(1-w0*f.subs(i, s+r))), sympy.Eq(s+r, N)),
    ((s*(s+2*r))/((2*(s+r)+1)*(2*(s+r+1)+1)), s >= N-r+1)
)

def get_gamma(l=l, v=v, n=N, r=r, n1=N1):
    if l == n1 + 1:
        return (
        (-1 + (sympy.exp(1)**(fi.subs({v: v})+fi_1.subs({v: v}))*1j/2) *
        ( (1+sympy.re(v)**2-sympy.im(v)**2)**2 + 4*sympy.re(v)**2*sympy.im(v)**2)**0.5  ) / 2 + 1
    )
    elif l.is_Number:
        gamma_l_plus_1 = get_gamma(l+1, v, n, r, n1)
    else:
        gamma_l_plus_1 = sympy.Function('gamma')(l+1, v, n, r, n1)

    return 1+a.subs( {l:l, v:v, r:r, w0:w0, n:n} ) / gamma_l_plus_1


def get_gamma_x(s=s, v=v, r=r, w0=w0):
    if s == N-r+1:
        return get_gamma(1, v, N, r)
    elif s.is_Number and r.is_Number:
        gamma_x_s_plus_1 = get_gamma_x(s+1, v, r, w0)
    else:
        gamma_x_s_plus_1 = sympy.Function('gamma_x')(s+1, v, r, w0)
    return 1 + q_x.subs( {s:s, r:r, w0:w0} )*v/gamma_x_s_plus_1


h_x = (
    1j*v*epsilon_x.subs( {s:alpha, r:r} ) *
    (kappa_x.subs( {s:alpha, r:r, w0:w0} )*get_gamma_x(alpha, v**2, r, w0))**(-1)
)

a = ((n+l)**2-r**2)*v**2/((2*(n+l)+1)*(2*(n+l+1)+1))


def get_Px(i=s+r+1, r=r, mu=mu):
    if i == r-1:
        return 0
    elif i == r:
        return (1-mu**2)**(r/2) * (2)**(-r) * sympy.sqrt(sympy.factorial(2*r)/sympy.factorial(r)**2)
    elif i.is_Number:
        Px_i_minus_1 = get_Px(i-1, r, mu)
        Px_i_minus_2 = get_Px(i-2, r, mu)
    else:
        Px_i_minus_1 = sympy.Function('Px')(i-1, r, mu)
        Px_i_minus_2 = sympy.Function('Px')(i-2, r, mu)

    return (
        ( 2*i*mu*Px_i_minus_1 - sympy.sqrt((i-r-1)*(i+r-1)) * Px_i_minus_2 ) *
        ( (i-r)*(i+r) )**(-1/2)
    )

psi_x = sympy.Piecewise(
    (1, sympy.Eq(j, 0)), 
    (sympy.Product(h_x.subs({l:l, v:v, r:r, w0:w0}), (l, 1, j)), True)
)

def get_nu_x(alpha=alpha, v=v, w0=w0, mu_1=mu_1, r=r):
    if alpha == N-r+1:
        return 0
    elif alpha == N-r:
        return (
            ((2*N+1)*f.subs({i:N})*get_Px(alpha+r, r, mu_1)) /
            (kappa_x.subs( {s:N-r, r:r, w0:w0} )*get_gamma_x(N-r, v**2, r, w0))
        )
    elif alpha.is_Number:    
        nu_x_alpha_plus_1 = get_nu_x(alpha+1, v, w0, mu_1)
    else:
        nu_x_alpha_plus_1 = sympy.Function('nu_x')(alpha+1, v, w0, mu_1)

    return (
        (kappa_x.subs( {s:alpha, r:r, w0:w0} )*get_gamma_x(alpha, v**2, r, w0))**(-1)*
        (1j*v*epsilon_x.subs( {s:alpha+1, r:r} )*nu_x_alpha_plus_1 + 
            (2*(alpha+r)+1)*f.subs({i:alpha+r})*get_Px(alpha+r, r, mu_1)
        )
    )


def get_b(alpha=alpha, v=v, w0=w0, mu_1=mu_1, r=r):
    if alpha == 0 and N-r >= 0:
        return (
            (kappa_x.subs({s:sympy.Integer('0'), r:r, w0:w0})*get_gamma_x(sympy.Integer('0'), v**2, r, w0))**(-1)*
            sympy.Sum(
                psi_x.subs({j:j, v:v, r:r, w0:w0})*(2*(j+r)+1)*f.subs({i:j+r})*get_Px(j+r, r, mu_1),
                (j, 0, N-r)
            )
        )
    elif alpha.is_Number and 1 <= alpha <= N-r and N-r > 1:
        b_alpha_minus_1 = get_b(alpha-1, v, w0, mu_1, r)
    else:
        b_alpha_minus_1 = sympy.Function('b')(alpha-1, v, w0, mu_1)
    
    return (
        h_x.subs({alpha:alpha, v:v, r:r, w0:w0})*b_alpha_minus_1 +
        get_nu_x(alpha, v, w0, mu_1)
    )

F = (0.5)*sympy.Sum(
    (2*(s+r+1))*get_b(s, v, w0, mu_1)*get_Px(s+r, r, mu),
    (s, 0, N)
)   
