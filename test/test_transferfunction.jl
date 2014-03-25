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
@test A + 1 == 1 + A == tf([1, 3], [1, 2])
@test A - 1 == tf([-1, -1], [1, 2])
@test 1 - A == tf([1, 1], [1, 2])
@test A*2 == 2*A == tf([2], [1, 2])
@test A/2 == tf([1], [2, 4])
@test 2/A == tf([2, 4], [1])

#Ones and Zeros
@test A - A == zero(A)
@test A*0 == zero(A)
@test A/A == one(A)
@test zero(A) == tf([0], [1])
@test one(A) == tf([1], [1])

#Indexing
@test A[1] == A.matrix[1]
@test size(A) == size(A.matrix)
@test length(A) == length(A.matrix) 

#Printing
@test string(A) == "  1.0\n-------\ns + 2.0"
