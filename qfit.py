#!/usr/bin/env python

import sys, numpy

A = open(sys.argv[1]).readlines()
p = len(A[2].split(".")[-2])+1
a,b,c = 20+p,20+2*p,20+3*p
X = numpy.array([(float(i[20:a]),float(i[a:b]),float(i[b:c])) for i in A[2:2+int(A[1])]])

B = open(sys.argv[2]).readlines()
p = len(B[2].split(".")[-2])+1
a,b,c = 20+p,20+2*p,20+3*p
Y = numpy.array([(float(i[20:a]),float(i[a:b]),float(i[b:c])) for i in B[2:2+int(B[1])]])

Y = Y[:X.shape[0],:]

m  = Y.mean(axis=0)
X -= X.mean(axis=0)

U,L,V = numpy.linalg.svd(numpy.dot(X.T,Y-m))

R = numpy.dot(U,V)
R[:,2] = numpy.cross(R[:,0],R[:,1])

print A[0],A[1],
for i,(x,y,z) in zip(A[2:],(numpy.dot(X,R)+m).tolist()):
    print i[:20]+("%8.3f%8.3f%8.3f"%(x,y,z))
print A[2+int(A[1])],
