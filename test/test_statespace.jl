using Base.Test
include("../src/Control.jl")
using Control

A = ss([1], [2], [3], [4])
B = ss([1 0; 6 1], [2 8]', [12 3], [16])
C = ss([-0.21 0.2; 0.2 -0.21], [0.01 0; 0 0.01], [1 0; 0 1], [0 0; 0 0])

#SS & SS Operations
@test A + B == ss([1 0 0; 0 1 0; 0 6 1], [2 2 8]', [3 12 3], [20])
@test A - B == ss([1 0 0; 0 1 0; 0 6 1], [2 2 8]', [3 -12 -3], [-12])
@test A * B == ss([1 0 0; 6 1 0; 24 6 1], [2 8 32]', [48 12 3], [64])
@test -A == ss([1], [2], [-3], [-4])

#SS & Constant Operations
ss([1], [2], [3], [4])
@test A + 1 == 1 + A == ss([1], [2], [3], [5])
@test A - 1 == ss([1], [2], [3], [3])
@test 1 - A == ss([1], [2], [-3], [-3])
@test A*2 == 2*A == ss([1], [2], [6], [8])
@test A/2 == ss([1], [2], [1.5], [2])

#Printing
@test sprint(show, B) == "StateSpace:\nA = \n1\t0\n6\t1\n\nB = \n2\n8\n\nC = \n12\t3\n\nD = \n16\n" 

#Test Errors
@test_throws A*C        #I/O dimension mismatch
@test_throws 1/A        #SS division not implemented
@test_throws A/B        #SS division not implemented
