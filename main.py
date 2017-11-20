#!/usr/bin/env python3
# -*-coding: utf-8 -*-

import sympy
import b_system


def main():
    sympy.init_printing()

    sympy.pprint(b_system.f)
    sympy.pprint(b_system.epsilon_x)
    sympy.pprint(b_system.kappa_x)
    sympy.pprint(b_system.fi)
    sympy.pprint(b_system.fi_1)
    sympy.pprint(b_system.q_x)
    sympy.pprint(b_system.get_gamma())
    sympy.pprint(b_system.get_gamma_x())
    sympy.pprint(b_system.h_x)
    sympy.pprint(b_system.a)
    sympy.pprint(b_system.psi_x)
    sympy.pprint(b_system.get_Px())
    sympy.pprint(b_system.get_nu_x())


if __name__ == '__main__':
    main()
