#####################################################################
##                      Data Type Declarations                     ##
#####################################################################

type StateSpace{T<:FloatingPoint} <: Sys
    A::Matrix{T}
    B::Matrix{T}
    C::Matrix{T}
    D::Matrix{T}
    states::Integer
    inputs::Integer
    outputs::Integer

    function StateSpace(A::Matrix{T}, B::Matrix{T}, C::Matrix{T}, D::Matrix{T})
        states = size(A)[1]
        inputs = size(B)[2]
        outputs = size(C)[1]

        if size(A)[2] != states
            error("A must be square")
        elseif size(B)[1] != states
            error("B must have the same row size as A")
        elseif size(C)[2] != states
            error("C must have the same column size as A")
        elseif inputs != size(D)[2]
            error("D must have the same column size as B")
        elseif outputs != size(D)[1]
            error("D must have the same row size as C")
        end

        new(A, B, C, D, states, inputs, outputs)
    end
end
StateSpace{T<:FloatingPoint}(A::Matrix{T}, B::Matrix{T}, C::Matrix{T}, D::Matrix{T}) = StateSpace{T}(A, B, C, D)

#####################################################################
##                      Constructor Functions                      ##
#####################################################################

## Fast Constructor for matrices of same type ##
ss{T<:FloatingPoint}(A::Matrix{T}, B::Matrix{T}, C::Matrix{T}, D::Matrix{T}) = StateSpace{T}(A, B, C, D)

## Constructor for 2x2 arrays of any real type ##
function ss{T1<:Real, T2<:Real, T3<:Real, T4<:Real}(A::Matrix{T1}, B::Matrix{T2}, C::Matrix{T3}, D::Matrix{T4})
    # Must be float, will always result in R = at least Float16
    R = promote_type(T1, T2, T3, T4, Float16)
    nA = convert(Matrix{R}, A)
    nB = convert(Matrix{R}, B)
    nC = convert(Matrix{R}, C)
    nD = convert(Matrix{R}, D)
    return StateSpace(nA, nB, nC, nD)
end

## General constructor, for 1D or 2D arrays of any real type ##
function ss{T1<:Real, T2<:Real, T3<:Real, T4<:Real}(A::Array{T1}, B::Array{T2}, C::Array{T3}, D::Array{T4})
    #Only check for 1x1 arrays. Anything else could be misinterpretted in
    #terms of transpose or not. Better to error than to assume.
    if ndims(A) == 1 && size(A)[1] == 1
        A = reshape(A, 1, 1)
    elseif ndims(A) > 2
        error("Matrices must not be larger than 2D")
    end
    if ndims(B) == 1 && size(B)[1] == 1
        B = reshape(B, 1, 1)
    elseif ndims(B) > 2
        error("Matrices must not be larger than 2D")
    end
    if ndims(C) == 1 && size(C)[1] == 1
        C = reshape(C, 1, 1)
    elseif ndims(C) > 2
        error("Matrices must not be larger than 2D")
    end
    if ndims(D) == 1 && size(D)[1] == 1
        D = reshape(D, 1, 1)
    elseif ndims(D) > 2
        error("Matrices must not be larger than 2D")
    end
    return ss(A, B, C, D)
end

#####################################################################
##                         Math Operators                          ##
#####################################################################

#TODO: Implement / for StateSpace/StateSpace
#----->Look at how OCTAVE implements this.

## EQUALITY ##
# Note: Due to multiple ss reps for the same i/o system, without
# trivial inequality does not necessarily mean the two systems aren't
# equivalent.
function ==(s1::StateSpace, s2::StateSpace)
    if s1.states != s2.states || s1.inputs != s2.inputs || s1.outputs != s2.outputs
        return false
    else
        return (s1.A == s2.A && s1.B == s2.B && s1.C == s2.C && s1.D == s2.D)
    end
end

## ADDITION ##
function +(s1::StateSpace, s2::StateSpace)
    #Ensure systems have same dimensions
    if s1.inputs != s2.inputs || s1.outputs != s2.outputs
        error("Systems have different shapes.")
    end

    A = [s1.A zeros(s1.states, s2.states); 
         zeros(s2.states, s1.states) s2.A]
    B = [s1.B ; s2.B]
    C = [s1.C s2.C]
    D = s1.D + s2.D

    return ss(A, B, C, D)
end

+(s1::StateSpace, s2::Real) = ss(s1.A, s1.B, s1.C, s1.D + s2)
+{T<:Real}(s2::T, s1::StateSpace) = +(s1, s2)

## SUBTRACTION ##
-(s1::StateSpace, s2::StateSpace) = +(s1, -s2)
-(s::StateSpace, n::Real) = ss(s.A, s.B, s.C, s.D - n)
-(n::Real, s::StateSpace) = +(-s, n)

## NEGATION ##
-(s1::StateSpace) = ss(s1.A, s1.B, -s1.C, -s1.D)

## MULTIPLICATION ##
function *(s1::StateSpace, s2::StateSpace)
    #Check dimension alignment
    if s1.inputs != s2.outputs
        error("A*B: # inputs A must equal # outputs B")
    end

    A = [s2.A  zeros(s2.states, s1.states); 
         s1.B*s2.C   s1.A]
    B = [s2.B ; s1.B*s2.D]
    C = [s1.D*s2.C   s1.C]
    D = s1.D * s2.D

    return ss(A, B, C, D)
end

*(s::StateSpace, n::Real) = ss(s.A, s.B, s.C*n, s.D*n)
*(n::Real, s::StateSpace) = *(s, n)

## DIVISION ##
function /(s1::StateSpace, s2::StateSpace)
    error("NotImplementedError: Division by StateSpace objects isn't
    implemented yet")
end

function /(n::Real, s::StateSpace)
    error("NotImplementedError: Division by StateSpace objects isn't
    implemented yet")
end

/(s::StateSpace, n::Real) = ss(s.A, s.B, s.C/n, s.D/n)

#####################################################################
##                        Display Functions                        ##
#####################################################################

function print(io::IO, s::StateSpace)
    print(io, "A = \n$(string(s.A))\n")
    print(io, "B = \n$(string(s.B))\n")
    print(io, "C = \n$(string(s.C))\n")
    print(io, "D = \n$(string(s.D))")
end

function show(io::IO, s::StateSpace)
    print(io, "StateSpace:\n")
    print(io, string(s))
end
