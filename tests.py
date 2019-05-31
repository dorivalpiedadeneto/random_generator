#! /usr/bin/python3
# This module contains tests for the random_generator module
# Tests:
# 01 - tests if intercepts works using floats to represent segments
# 02 - tests if intercepts works using arrays to represent segments 
from generator import intercepts, points_inside_boundary
from numpy import array as nparray

def test01():
    # 1sd segment (boundary)
    s1 = ((0.0,0.0),(0.0,1.0))
    # segments
    s11 = ((-1.0,0.0),(1.0,0.0))
    s12 = ((-1.0,0.5),(1.0,0.5))
    s13 = ((-1.0,1.0),(1.0,1.0))
    s14 = ((-1.0,-1.0e-6),(1.0,-1.0e-6))
    s15 = ((-1.0,1.000001),(1.0,1.00001))
    s_ = (s11,s12,s13,s14,s15)
    for s in s_:
        print(' - Testing interception between {} and {}'.format(s,s1))
        if intercepts(s,s1):
            print('   Intercepts!')
        else:
            print('    Does not intercept!')

def test02():
    s1 = ((0.0,0.0),(0.0,1.0))
    # segments
    xi = nparray((-1.0,-1.0,-1.0,-1.0,-1.0))
    xf = nparray((1.0,1.0,1.0,1.0,1.0))
    yi = nparray((0.0,0.5,1.0,-1.0e-6,1.00001))
    yf = nparray((0.0,0.5,1.0,-1.0e-6,1.00001))
    s_ = ((xi,xf),(yi,yf))
    print('Testing interception of {} and segments from array {}'.format(s1,s_))
    print(intercepts(s_,s1))

def test03():
    pts_ = ((0.0,0.0),(2.0,2.0))
    bvs_ = ((-1.0,-1.0),
            ( 1.0,-1.0),
            ( 1.0, 1.0),
            (-1.0, 1.0),
            (-1.0,-1.0))
    res = points_inside_boundary(pts_,bvs_)
    print('Points inside area:')
    print(res)

if __name__ == "__main__":

    test01()
    test02()
    test03()
