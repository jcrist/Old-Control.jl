## Module for dealing with State Space system representation ##

#Types:
#--> StateSpace

#Constructors:
#--> ss(A, B, C, D)

#Methods:
#--> show
#--> string(TransferFunction)
#--> Operators: +, -, *, /

module statespace

export ss, StateSpace

import Base: length, getindex, show, string, print
import Base: *, /, +, -

#####################################################################
##                      Data Type Declarations                     ##
#####################################################################

type StateSpace
    A::Array{Float64, 2}
    B::Array{Float64, 2}
    C::Array{Float64, 2}
    D::Array{Float64, 2}
    states::Integer
    inputs::Integer
    outputs::Integer
    function StateSpace(A::Array{Float64, 2}, B::Array{Float64, 2}, C::Array{Float64, 2}, D::Array{Float64, 2})
        states = size(A)[1]
        inputs = size(B)[2]
        outputs = size(C)[1]

        #Perform checks on input, ensure shapes align
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

        #Copy matrices to prevent change by reference behavior
        A = copy(A)
        B = copy(B)
        C = copy(C)
        D = copy(D)

        #Return a StateSpace Type
        new(A, B, C, D, states, inputs, outputs)
    end
end

#####################################################################
##                      Constructor Functions                      ##
#####################################################################

#TODO: There has to be a better way to design these constructors!!?!

function ss{T1<:Real, T2<:Real, T3<:Real, T4<:Real}(A::Array{T1}, B::Array{T2}, C::Array{T3}, D::Array{T4})
    ## General constructor, for any array of any real type ##

    #Only check for 1x1 arrays. Anything else could be misinterpretted in
    #terms of transpose or not. Better to error than to assume.
    if ndims(A) == 1 && size(A)[1] == 1
        A = reshape(A, 1, 1)
    elseif ndims(A) > 2
        error("Matrices must be not be larger than 2D")
    end
    if ndims(B) == 1 && size(B)[1] == 1
        B = reshape(B, 1, 1)
    elseif ndims(B) > 2
        error("Matrices must be not be larger than 2D")
    end
    if ndims(C) == 1 && size(C)[1] == 1
        C = reshape(C, 1, 1)
    elseif ndims(C) > 2
        error("Matrices must be not be larger than 2D")
    end
    if ndims(D) == 1 && size(D)[1] == 1
        D = reshape(D, 1, 1)
    elseif ndims(D) > 2
        error("Matrices must be not be larger than 2D")
    end
    nA = convert(Array{Float64, 2}, A)
    nB = convert(Array{Float64, 2}, B)
    nC = convert(Array{Float64, 2}, C)
    nD = convert(Array{Float64, 2}, D)
    return StateSpace(nA, nB, nC, nD)
end

function ss{T1<:Real, T2<:Real, T3<:Real, T4<:Real}(A::Array{T1, 2}, B::Array{T2, 2}, C::Array{T3, 2}, D::Array{T4, 2})
    #General constructor, for 2x2 arrays of any real type
    nA = convert(Array{Float64, 2}, A)
    nB = convert(Array{Float64, 2}, B)
    nC = convert(Array{Float64, 2}, C)
    nD = convert(Array{Float64, 2}, D)
    return StateSpace(nA, nB, nC, nD)
end

#Fast constructor for arrays that are all Float64
ss(A::Array{Float64, 2}, B::Array{Float64, 2}, C::Array{Float64, 2}, D::Array{Float64, 2}) = StateSpace(A, B, C, D)


#####################################################################
##                         Math Operators                          ##
#####################################################################

#TODO: Implement / for StateSpace/StateSpace
#----->Look at how OCTAVE implements this.
#TODO: Some of these could be more efficient, many use intermediate
#----->terms, and pass on to an already defined function. As speed for
#----->these operations isn't imperative, this is ok, for now.

## Addition ##
function +(self::StateSpace, other::StateSpace)
    ## Add two state space objects ##

    #Ensure systems have same dimensions
    if self.inputs != other.inputs || self.outputs != other.outputs
        error("Systems have different shapes.")
    end

    A = [self.A zeros(size(self.A)); zeros(size(other.A)) other.A]
    B = [self.B ; other.B]
    C = [self.C other.C]
    D = self.D + other.D

    return StateSpace(A, B, C, D)
end

function +{T<:Real}(self::StateSpace, other::T)
    ## Add a number to a state space object ##
    return StateSpace(self.A, self.B, self.C, self.D + other)
end

+{T<:Real}(other::T, self::StateSpace) = +(self, other)

## Subtraction ##
function -(self::StateSpace, other::StateSpace)
    ## Subtract two state space objects ##
    return +(self, -other)
end

function -{T<:Real}(self::StateSpace, other::T)
    ## Subtract a number from a state space object ##
    return StateSpace(self.A, self.B, self.C, self.D - other)
end

-{T<:Real}(other::T, self::StateSpace) = +(-self, other)


## Negation ##
function -(self::StateSpace)
    return StateSpace(self.A, self.B, -self.C, -self.D)
end

## Multiplication ##
function *(self::StateSpace, other::StateSpace)
    ## Multiply two state space objects ##

    #Check dimension alignment
    if self.inputs != other.outputs
        error("A*B: # inputs A must equal # outputs B")
    end

    #Multiply/Concatenate matrices to form resulting system
    A = [[other.A   zeros(size(other.A))] ; [self.B*self.C   self.A]]

    B = [other.B ; self.B*other.D]

    C = [self.D*other.C   self.C]

    D = self.D * other.D

    return StateSpace(A, B, C, D)
end

function *{T<:Real}(self::StateSpace, other::T)
    ## Multiply a state space object and a scalar ##
    return StateSpace(self.A, self.B, self.C * other, self.D * other)
end

*{T<:Real}(other::T, self::StateSpace) = *(self, other)

## Division ##
function /(self::StateSpace, other::StateSpace)
    ## Divide two state space objects ##

    error("NotImplementedError: Division by StateSpace objects isn't
    implemented yet")
end

function /{T<:Real}(self::StateSpace, other::T)
    ## Divide a state space object with a scalar ##
    return StateSpace(self.A, self.B, self.C / other, self.D / other)
end

function /{T<:Real}(other::T, self::StateSpace)
    ## Divide a scalar with a state space object ##
    error("NotImplementedError: Division by StateSpace objects isn't
    implemented yet")
end


#####################################################################
##                        Display Functions                        ##
#####################################################################

#TODO: String representation of matrices could be prettier. Currently
#      relying on standard array string function, which lacks braces.
#      Possible format would use | on edges.

function string(self::StateSpace)
    ## String representation of a state space system ##    

    str = "A = \n$(string(self.A))"
    str = "$str\nB = \n$(string(self.B))"
    str = "$str\nC = \n$(string(self.C))"
    str = "$str\nD = \n$(string(self.D))\n"

    return str
end

function show(io::IO, self::StateSpace)
    print("StateSpace:\n")
    print(string(self))
end

end     #Module
