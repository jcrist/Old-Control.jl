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
@test 1/A == ss([-0.5], [0.5], [-0.75], [0.25])
@test 1/(1/A) == A
@test A/B == ss([1 0 0; 0.375 -0.5 -0.375; 1.5 0 -0.5], [2 0.5 2]', [0.1875 -0.75 -0.1875], [0.25])

#Printing
@test sprint(show, C) == "StateSpace:\nA = \n[-0.21 0.2\n 0.2 -0.21]\nB = \n[0.01 0.0\n 0.0 0.01]\nC = \n[1.0 0.0\n 0.0 1.0]\nD = \n[0.0 0.0\n 0.0 0.0]"

#Test Errors
@test_throws ErrorException A*C                                        #I/O dimension mismatch
@test_throws ErrorException ss([1 0; 6 1], [2; 8], [12 3], [16])       #Constructor with 1D array with more than 1 entry
@test_throws ErrorException 1/C                                        #Non-invertible D matrix
