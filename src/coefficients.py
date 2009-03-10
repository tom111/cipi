#!/usr/bin/python

# This program will find the coefficients in an energy expansion of a
# probability vector

import sys # Mainly for file IO
from math import *
from array import * #For fast arrays

def readPVector (filename):
    """ reads a probability vector from the file given in filename"""
    try:
        f = open(filename,'r')
    except IOError:
        print 'File could not be opened'
        exit(1);
    lines = f.readlines()
    splitlines = []
    for line in lines:
        splitlines += [line.split()]

    return array ('d', [float(x[1]) for x in splitlines])

def bitString(s):
    return str(s) if s<=1 else bitString(s>>1) + str(s&1)

def parity (s):
    """ Returns the parity of a bitstring"""
    su = 0;
    for b in s:
        if b=='1':
            su=su+1
    return su%2

def bitWeight (n):
    """ Returns the number of ones in a bitstring"""
    su = 0;
    for b in bitString(n):
        if b=='1':
            su=su+1
    return su;

def E(A,x):
    """ Returns the E(A,x), A,x given as integers """
    return (-1)**( bitWeight (A&x) )

def project(A,H):
    count=0
    su=0
    for h in H:
        su += E(A,count) *h
        count+=1
    return su

def bitStringtoSet (s):
    s = s[::-1] # Reverses the string
    se = []
    count = 0
    for b in s:
        count += 1
        if b=='1':
            se.append(count);
    return se
        

######## Main program starts here ##############

try:
    Pfile = sys.argv[1]
except IndexError:
    print 'Please give filename as first argument';
    exit(1);

P = readPVector(Pfile)

for p in P:
    if p<=0 or p>=1:
        print 'Error, support of P is not full'
        exit(1)

# Take the logarithm of the probabilities
P = [log (p) for p in P]
print P

# Iterating over the number of elements in P
n2 = len(P)

result = []
for a in range (n2):
    result.append([a, project(a,P)])

# Sorting the coefficients by bitweight
result.sort(key=lambda t:bitWeight(t[0]))

for c in result:
    print  ' : '.join( [str(bitStringtoSet(bitString(c[0]))),str( c[1])])
