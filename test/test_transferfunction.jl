using Base.Test
include("../src/Control.jl")
using Control

A = tf([1], [1, 2])
B = tf([1], [1, 2, 3])

#TF & TF Operations
@test A+B == tf([1, 3, 5], [1, 4, 7, 6])
@test A-B == tf([1, 1, 1], [1, 4, 7, 6])
@test A*B == tf([1], [1, 4, 7, 6])
@test A/B == tf([1, 2, 3], [1, 2])

#TF & Constant Operations
@test A + 1 == tf([1, 3], [1, 2])
@test A - 1 == tf([-1, -1], [1, 2])
@test A*2 == tf([2], [1, 2])
@test A/2 == tf([1], [2, 4])

#Printing
@test string(A) == "  1\n-----\ns + 2\n\n" 
