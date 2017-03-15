#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Li Zhu < zhulipresent@gmail.com >


import os
import sys
import math

def writekp(kgrid):
    """
    
    Arguments:
    - `kgrid` : Kmesh
    """
    # read the lattice

    def dot(x, y):
        return x[0]*y[0] + x[1]*y[1] + x[2]*y[2]

    def cross(x, y):
        z1 = x[1]*y[2] - x[2]*y[1]
        z2 = x[2]*y[0] - x[0]*y[2]
        z3 = x[0]*y[1] - x[1]*y[0]
        return [z1, z2, z3]
    def kmf(kgrid, gi):
        kd = int(gi/kgrid/2.0/math.pi)
        if kd == 0: kd = 1
        dd = gi/kd/2.0/math.pi
        if dd >= kgrid:
            for i in range(0, 10):
                kd += i
                dd = gi/kd/2.0/math.pi
                if dd <= kgrid: break
        return kd 
                    
    f = open('POSCAR')
    pp = []
    try:
        for line in f:
            pp.append(line.split())
    finally:
        f.close()
    l = []
    for item in pp[2:5]:
        l.append(map(float, item))
    # real Lattice Parameters
    ra = math.sqrt(l[0][0]**2 + l[0][1]**2 + l[0][2]**2)
    rb = math.sqrt(l[1][0]**2 + l[1][1]**2 + l[1][2]**2)
    rc = math.sqrt(l[2][0]**2 + l[2][1]**2 + l[2][2]**2)

    c = cross(l[1], l[2])
    volume = dot(l[0], c)
    g = []
    g1 = [ 2.0 * math.pi * item / volume for item in c]
    c = cross(l[2], l[0])
    g2 = [ 2.0 * math.pi * item / volume for item in c]
    c = cross(l[0], l[1])
    g3 = [ 2.0 * math.pi * item / volume for item in c]
    g = [g1, g2, g3]
    
    rl = []
    for i in range(0, 3):
        rl.append(math.sqrt(dot(g[i],g[i])))
    kmesh = []
    for i in range(0, 3):
        kmesh.append(kmf(kgrid, rl[i]))
    f = open('KPOINTS', 'w')
    f.write('A\n0\nG\n')
    f.write('%2d %2d %2d\n' % tuple(kmesh))
    f.write('%2d %2d %2d\n' % (0,0,0))
    # debug
    #print 'volume: ', volume
    #print 'rl ', rl
    #print 'kp ', kmesh
    # end debug

if __name__ == '__main__':
    args = sys.argv
    try:
        kgrid = float(args[1])
    except:
        kgrid = 0.1
    writekp(kgrid)
