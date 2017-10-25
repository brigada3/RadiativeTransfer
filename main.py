#!/usr/bin/env python3
# -*-coding: utf-8 -*-

from b_system import *


def main():
    print("epsilon_x:", epsilon_x(1000, 1000))
    print("kappa_x:", kappa_x(1000, 1000, 1))
    print("fi:", fi(2-3j))
    print("fi_1:", fi_1(2-3j))
    print("q_x:", q_x(1000, 0, 0))
    print("h_x:", h_x(2000, 2j+3, 1, 0))


if __name__ == '__main__':
    main()
    