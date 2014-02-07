#Module for working with Transfer Functions

#Types:
#--> TransferFunction

#Constructors:
#--> tf(num, den)

#Methods:
#--> show
#--> string(TransferFunction)
#--> Operators: +, -, *, /

module XferFcn

export tf, TransferFunction

import Base: length, getindex, show, string, print
import Base: *, /, +, -
import Polylib: polymul, polydiv, polyadd, polytostring, polyfracsimp, trimzeros


#TODO: Where should these go?
#Used for determining if value is ~0
eps{T}(::Type{T}) = convert(T,0)
eps{F<:FloatingPoint}(x::Type{F}) = Base.eps(F)


#####################################################################
##                      Data Type Declarations                     ##
#####################################################################

type TransferFunction
    #Datatype for SISO and MIMO transfer functions
    num::Array{Vector{Float64}, 2}
    den::Array{Vector{Float64}, 2}
    inputs::Integer
    outputs::Integer

    #Inner constructor for input validation, and determining inputs/outputs
    function TransferFunction(num::Array{Vector{Float64}, 2}, den::Array{Vector{Float64}, 2})
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

#TODO: Allow for num, den of any type, convert to float

function tf(num::Vector{Float64}, den::Vector{Float64})
    #Create SISO system
    narr = Array(Vector{Float64}, 1, 1)
    narr[1,1] = num
    darr = Array(Vector{Float64}, 1, 1)
    darr[1,1] = den
    TransferFunction(narr, darr)
end


#####################################################################
##                         Math Operators                          ##
#####################################################################

#TODO: Implement *, /
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
    ## Add a transfer function to a number ##
    return +(self, tf([other], [1.0]))
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
    return -(self, tf([other], [1.0]))
end

-{T<:Real}(other::T, self::TransferFunction) = -(self, other)

## NEGATION ##
function -(self::TransferFunction)
    ## Negate a transfer function ##
    #Preallocate space for num, deepcopy den
    num = Array(Vector{Float64}, self.outputs, self.inputs)
    den = deepcopy(self.den)
    for o=1:self.outputs
        for i=1:self.inputs
            num[o,i] = -self.num[o,i]
        end
    end
    return TransferFunction(num, den)
end

## HELPERS ##
function _addSISO(num1::Vector{Float64}, den1::Vector{Float64}, num2::Vector{Float64}, den2::Vector{Float64})
    #Helper function for adding 2 SISO systems together
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
    print("TransferFunction\n")
    print(string(tf))
end

end     #Module
