#Module for working with Transfer Functions

#Types:
#--> TransferFunction

#Constructors:
#--> tf(num, den)

#Conversions:
#--> ss2tf(sys)

#Methods:
#--> show
#--> string(TransferFunction)
#--> Operators: +, -, *, /

module transferfunction

export tf, TransferFunction, ss2tf

import Base: length, getindex, show, string, print
import Base: *, /, +, -
import polylib: polymul, polydiv, polyadd, polytostring, polyfracsimp, 
                trimzeros
import statespace: StateSpace

import Control: Sys
import Slicot.Simple: tb04ad


#####################################################################
##                      Data Type Declarations                     ##
#####################################################################

type TransferFunction <: Sys
    ## Datatype for SISO and MIMO transfer functions ##
    num::Array{Vector{Float64}, 2}
    den::Array{Vector{Float64}, 2}
    inputs::Integer
    outputs::Integer

    function TransferFunction(num::Array{Vector{Float64}, 2}, den::Array{Vector{Float64}, 2})
        ## Inner constructor for input validation, and determining inputs/outputs ##
        #Validate input and output dimensions match
        out_n, in_n = size(num)
        out_d, in_d = size(den)
        if out_n != out_d
            error("num and den output dimensions must match")
        end
        outputs = out_n

        if in_n != in_d
            error("num and den input dimensions must match")
        end
        inputs = in_n

        ## Remove leading zeros on num and den ##
        #Data is deepcopied to prevent mutating the calling arrays
        #TODO: Determine if deepcopy needed. Should be able to allocate space,
        #and just return trimmed arrays into location.
        new_den = deepcopy(den)
        new_num = deepcopy(num)
        for o=1:outputs
            for i=1:inputs
                den_temp = trimzeros(new_den[o,i])
                if den_temp == [0.0]
                    #This is a denominator with zero value, throw error
                    error("Input $i, output $o has a zero denominator")
                end
                num_temp = trimzeros(new_num[o,i])
                if num_temp == [0.0]
                    #The numerator is zero, make the denominator 1
                    #assignment to the array occurs here, as no need to simplify
                    new_num[o,i] = num_temp
                    new_den[o,i] = [1.0]
                    #No need to simplify the fraction, continue to next input
                    continue
                end
                #Simplify the resulting fraction
                new_num[o,i], new_den[o,i] = polyfracsimp(num_temp, den_temp)
            end
        end
        new(new_num, new_den, inputs, outputs)
    end
end


#####################################################################
##                      Constructor Functions                      ##
#####################################################################

function tf{T<:Real, S<:Real}(num::Vector{T}, den::Vector{S})
    ## Create SISO system ##
    narr = Array(Vector{Float64}, 1, 1)
    narr[1,1] = convert(Vector{Float64}, num)
    darr = Array(Vector{Float64}, 1, 1)
    darr[1,1] = convert(Vector{Float64}, den)
    TransferFunction(narr, darr)
end

function tf(num::Vector{Float64}, den::Vector{Float64})
    ## Create SISO system ##
    #Version for already Float64, avoids type conversion
    narr = Array(Vector{Float64}, 1, 1)
    narr[1,1] = num
    darr = Array(Vector{Float64}, 1, 1)
    darr[1,1] = den
    TransferFunction(narr, darr)
end

#####################################################################
##                      Conversion Functions                       ##
#####################################################################
function ss2tf(sys::StateSpace)
    ## Convert a StateSpace to a TransferFunction ##
    tfout = tb04ad('R', sys.states, sys.inputs, sys.outputs, sys.A,
                    sys.B, sys.C, sys.D)

    #Allocate space for the num and den arrays
    num = Array(Vector{Float64}, sys.outputs, sys.inputs)
    den = Array(Vector{Float64}, sys.outputs, sys.inputs)

    for o=1:sys.outputs
        for i=1:sys.inputs
            n1, n2, n3 = size(tfout[7][o,i,:])
            num[o,i] = reshape(tfout[7][o, i, :], n3)
            m,n = size(tfout[6][o, :])
            den[o,i] = reshape(tfout[6][o, :], n)
        end
    end
    return TransferFunction(num, den)
end

function ss2tf(A::Array{Float64,2}, B::Array{Float64,2}, C::Array{Float64,2}, D::Array{Float64,2})
    ## Convert system matrices into a TransferFunction ##

    #Calculate necessary dimensions
    #Note that no size verification occurs here. As tb04ad does this already,
    #it would be redundant to do it again. Perhaps this may be changed
    #if needed.
    states = size(A)[1]
    inputs = size(B)[2]
    outputs = size(C)[1]
    tfout = tb04ad('R', states, inputs, outputs, A, B, C, D)

    #Allocate space for the num and den arrays
    num = Array(Vector{Float64}, sys.outputs, sys.inputs)
    den = Array(Vector{Float64}, sys.outputs, sys.inputs)

    for o=1:sys.outputs
        for i=1:sys.inputs
            n1, n2, n3 = size(tfout[7][o,i,:])
            num[o,i] = reshape(tfout[7][o, i, :], n3)
            m,n = size(tfout[6][o, :])
            den[o,i] = reshape(tfout[6][o, :], n)
        end
    end
    return TransferFunction(num, den)
end

#####################################################################
##                         Math Operators                          ##
#####################################################################

#TODO: Refactor for less code, using parametric types.
#----> Should be able to call operator, and have operator attempt
#----> conversion to TF if other is not a TF.


## ADDITION ##
function +(self::TransferFunction, other::TransferFunction)
    ## Add two transfer function objects (parallel connection) ##

    #Check that the input-output sizes are consistent
    if self.inputs != other.inputs
        error("The first summand has %i input(s), but the second has %i.", self.inputs, other.inputs)
    end
    if self.outputs != other.outputs
        error("The first summand has %i output(s), but the second has %i.", self.outputs, other.outputs)
    end

    #Preallocate the numerator and denominator of the sum
    num = Array(Vector{Float64}, self.outputs, self.inputs)
    den = Array(Vector{Float64}, self.outputs, self.inputs)

    #Perform SISO addition for each input & output
    for o=1:self.outputs
        for i=1:self.inputs
            num[o,i], den[o, i] = _addSISO(self.num[o,i], self.den[o,i], other.num[o,i], other.den[o,i])
        end
    end

    return TransferFunction(num, den)
end

function +{T<:Real}(self::TransferFunction, other::T)
    ## Add a number to a transfer function ##
    temp = tf([other], [1.0])
    #Ensure input-output match
    #this is kinda hacky, but it works for the checks we need
    temp.inputs = self.inputs
    temp.outputs = self.outputs
    return +(self, temp)
end

+{T<:Real}(other::T, self::TransferFunction) = +(self, other)
    
## SUBTRACTION ##
function -(self::TransferFunction, other::TransferFunction)
    ## Subtract two transfer function objects (parallel connection) ##

    #Check that the input-output sizes are consistent
    if self.inputs != other.inputs
        error("The first summand has %i input(s), but the second has %i.", self.inputs, other.inputs)
    end
    if self.outputs != other.outputs
        error("The first summand has %i output(s), but the second has %i.", self.outputs, other.outputs)
    end

    #Preallocate the numerator and denominator of the sum
    num = Array(Vector{Float64}, self.outputs, self.inputs)
    den = Array(Vector{Float64}, self.outputs, self.inputs)

    #Perform SISO addition for each input & output
    for o=1:self.outputs
        for i=1:self.inputs
            num[o,i], den[o, i] = _addSISO(self.num[o,i], self.den[o,i], -other.num[o,i], other.den[o,i])
        end
    end

    return TransferFunction(num, den)
end

function -{T<:Real}(self::TransferFunction, other::T)
    ## Subtract a number from a transfer function ##
    temp = tf([other], [1.0])
    #Ensure input-output match
    #this is kinda hacky, but it works for the checks we need
    temp.inputs = self.inputs
    temp.outputs = self.outputs
    return -(self, temp)
end

-{T<:Real}(other::T, self::TransferFunction) = +(-self, other)

## NEGATION ##
function -(self::TransferFunction)
    ## Negate a transfer function ##
    #Preallocate space for num
    num = Array(Vector{Float64}, self.outputs, self.inputs)
    for o=1:self.outputs
        for i=1:self.inputs
            num[o,i] = -self.num[o,i]
        end
    end
    return TransferFunction(num, self.den)
end

## MULTIPLICATION ##
function *(self::TransferFunction, other::TransferFunction)
    ## Multiply two transfer functions together ##

    #Check that the input-output sizes are consistent
    if self.inputs != other.outputs
        error(@sprintf("Input->Output Mismatch: C = A*B: A has %i inputs, B has %i outputs", self.inputs, other.outputs))
    end

    inputs = other.inputs
    outputs = self.outputs

    #Preallocate the numerator and denominator of the product
    num = Array(Vector{Float64}, self.outputs, self.inputs)
    den = Array(Vector{Float64}, self.outputs, self.inputs)

    #Temporary storage for the summands needed to find the (o, i)th element
    #of the product.
    num_summand = Array(Vector{Float64}, self.inputs)
    den_summand = Array(Vector{Float64}, self.inputs)

    for o=1:outputs
        for i=1:inputs
            for k=1:self.inputs
                num_summand[k] = polymul(self.num[o,k], other.num[k,i])
                den_summand[k] = polymul(self.den[o,k], other.den[k,i])
                num[o,i], den[o,i] = _addSISO([0.0], [1.0], num_summand[k], den_summand[k])
            end
        end
    end

    return TransferFunction(num, den)
end

function *{T<:Real}(self::TransferFunction, other::T)
    ## Multiply a number to a transfer function ##
    temp = tf([other], [1.0])
    #Ensure input-output match
    #this is kinda hacky, but it works for the checks we need
    temp.inputs = self.inputs
    temp.outputs = self.outputs
    return *(self, temp)
end

*{T<:Real}(other::T, self::TransferFunction) = *(self, other)

## DIVISION ##
function /(self::TransferFunction, other::TransferFunction)
    ## Divide two transfer functions ##
    #TODO: Implement division for MIMO systems

    if self.inputs > 1 || self.outputs > 1 || other.inputs > 1 || other.outputs > 1
        error("NotImplementedError: Division is currently only implemented for SISO systems")
    end

    num = polymul(self.num[1,1], other.den[1,1])
    den = polymul(self.den[1,1], other.num[1,1])

    return tf(num, den)
end

/{T<:Real}(self::TransferFunction, other::T) = /(self, tf([other], [1.0]))
/{T<:Real}(other::T, self::TransferFunction) = /(tf([other], [1.0]), self)

## HELPERS ##
function _addSISO(num1::Vector{Float64}, den1::Vector{Float64}, num2::Vector{Float64}, den2::Vector{Float64})
    ## Helper function for adding 2 SISO systems together ##
    num = polyadd(polymul(num1, den2), polymul(num2, den1))
    den = polymul(den1, den2)
    return num, den
end


#####################################################################
##                        Display Functions                        ##
#####################################################################

function string(tf::TransferFunction)
    ## String representation of the transfer function ##

    mimo = tf.inputs > 1 || tf.outputs > 1
    outstr = ""

    for i=1:tf.inputs
        for j=1:tf.outputs
            if mimo
                outstr = "$(outstr)Input $i to output $j\n"
            end

            #Convert the numerator and denominator to strings
            numstr = polytostring(tf.num[j,i])
            denstr = polytostring(tf.den[j,i])

            #Figure out the length of the separating line
            dashcount = max(length(numstr), length(denstr))

            #Center the numerator or denominator
            if length(numstr) < dashcount
                numstr = "$(repeat(" ", int((dashcount - length(numstr))/2)))$numstr"
            else
                denstr = "$(repeat(" ", int((dashcount - length(denstr))/2)))$denstr"
            end
            #Cat string with numerator, dashes, and denominator strings
            outstr = "$outstr$numstr\n$(repeat("-", dashcount))\n$denstr\n\n"
        end
    end

    return outstr
end

function show(io::IO, tf::TransferFunction)
    ## Print a transfer function representation in the REPL ##
    print("TransferFunction:\n")
    print(string(tf))
end

end     #Module
