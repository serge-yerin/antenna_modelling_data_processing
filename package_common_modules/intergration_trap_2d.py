'''
'''
import numpy as np

def intergration_trap_2d (a, b, c, d, m, n, Points):
    pi = 3.141593

    # a - the start integration limit by the first axis
    # b - the end integration limit by the first axis
    # c - the start integration limit by the second axis
    # d - the end integration limit by the second axis
    # n - number of points on 1 axis
    # m - number of points on 2 axis

    #Sum = 0.0
    #for i in range (m):
    #    for k in range (n):
    #        Sum = Sum + Points[i, k]

    Sum = np.sum(Points)

    IntegrResult = ((((b)-(a))*pi/180.0)/m) * ((((d)-(c))*pi/180.0)/n) * ((Points[a,c] + Points[b,d]) + Sum)
    return IntegrResult
