#####################################################################
##                      Data Type Declarations                     ##
#####################################################################

## User should just use TransferFunction
immutable SISOTransferFunction{T<:FloatingPoint} <: Sys
    ## Type for SISO Transfer Functions only. ##
    num::Poly{T}
    den::Poly{T}
    function SISOTransferFunction(num::Poly{T}, den::Poly{T})
        fact = gcd(num, den)
        num_s = num/fact
        den_s = den/fact
        #To stop the num and den from getting smaller each time, calculate
        #a scaling factor so that den[1] = 1
        scale = 1.0/den_s[1]
        num_s = num_s*scale
        den_s = den_s*scale
        return new(num_s, den_s)
    end
end
SISOTransferFunction{T<:FloatingPoint}(num::Poly{T}, den::Poly{T}) = SISOTransferFunction{T}(num, den)

type TransferFunction{T<:FloatingPoint} <: Sys
    ## Datatype for SISO and MIMO transfer functions ##
    matrix::Array{SISOTransferFunction{T}, 2}
    inputs::Integer
    outputs::Integer
end

function TransferFunction{T<:FloatingPoint}(num::Array{Vector{T}, 2}, den::Array{Vector{T}, 2})
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

    matrix = Array(SISOTransferFunction{T}, outputs, inputs)
    for o=1:outputs
        for i=1:inputs
            den_temp = Poly(den[o,i], 's')
            if den_temp == zero(den_temp)
                #This is a denominator with zero value, throw error
                error("Input $i, output $o has a zero denominator")
            end
            num_temp = Poly(num[o,i], 's')
            if num_temp == zero(num_temp) 
                #The numerator is zero, make the denominator 1
                den_temp = one(den_temp)
            end
            matrix[o,i] = SISOTransferFunction(num_temp, den_temp)
        end
    end
    TransferFunction(matrix, inputs, outputs)
end


#####################################################################
##                      Constructor Functions                      ##
#####################################################################

function tf{T<:Real, S<:Real}(num::Vector{T}, den::Vector{S})
    ## Create SISO system ##
    R = promote_type(promote_type(T, S), Float16)
    narr = Array(Vector{R}, 1, 1)
    narr[1,1] = convert(Vector{R}, num)
    darr = Array(Vector{R}, 1, 1)
    darr[1,1] = convert(Vector{R}, den)
    TransferFunction(narr, darr)
end

#####################################################################
##                      Conversion Functions                       ##
#####################################################################
promote_rule{T, S}(::Type{SISOTransferFunction{T}}, ::Type{SISOTransferFunction{S}}) = SISOTransferFunction{promote_type(T, S)}
promote_rule{T, S<:Real}(::Type{SISOTransferFunction{T}}, ::Type{S}) = SISOTransferFunction{promote_type(T, S)}
promote_rule{T, S}(::Type{TransferFunction{T}}, ::Type{TransferFunction{S}}) = TransferFunction{promote_type(T, S)}

convert{T}(::Type{SISOTransferFunction{T}}, t::SISOTransferFunction) = 
    SISOTransferFunction(convert(Poly{T}, t.num), convert(Poly{T}, t.den))
convert{T}(::Type{SISOTransferFunction{T}}, o::Real) = o*one(SISOTransferFunction{T})

if SLICOT_DEFINED
function ss2tf(ss::StateSpace)
    ## Convert a StateSpace to a TransferFunction ##
    tfout = tb04ad('R', ss.states, ss.inputs, ss.outputs, ss.A,
                    ss.B, ss.C, ss.D)

    #Allocate space for the num and den arrays
    num = Array(Vector{Float64}, ss.outputs, ss.inputs)
    den = Array(Vector{Float64}, ss.outputs, ss.inputs)

    for o=1:ss.outputs
        for i=1:ss.inputs
            n1, n2, n3 = size(tfout[7][o,i,:])
            num[o,i] = reshape(tfout[7][o, i, :], n3)
            m,n = size(tfout[6][o, :])
            den[o,i] = reshape(tfout[6][o, :], n)
        end
    end
    return TransferFunction(num, den)
end

function ss2tf(A::Array{Float64,2}, B::Array{Float64,2}, C::Array{Float64,2}, D::Array{Float64,2})
    ## Convert StateSpace to TransferFunction ##

    #Calculate necessary dimensions
    states = size(A)[1]
    inputs = size(B)[2]
    outputs = size(C)[1]
    tfout = tb04ad('R', states, inputs, outputs, A, B, C, D)

    #Allocate space for the num and den arrays
    num = Array(Vector{Float64}, ss.outputs, ss.inputs)
    den = Array(Vector{Float64}, ss.outputs, ss.inputs)

    for o=1:ss.outputs
        for i=1:ss.inputs
            n1, n2, n3 = size(tfout[7][o,i,:])
            num[o,i] = reshape(tfout[7][o, i, :], n3)
            m,n = size(tfout[6][o, :])
            den[o,i] = reshape(tfout[6][o, :], n)
        end
    end
    return TransferFunction(num, den)
end
end #SLICOT_DEFINED

#####################################################################
##                          Misc. Functions                        ##
#####################################################################

## ONE ##
one(t::SISOTransferFunction) = SISOTransferFunction(one(t.num), one(t.num))
one{T}(::Type{SISOTransferFunction{T}}) = SISOTransferFunction(Poly([one(T)], :s), Poly([one(T)], :s))
function one{T}(t::TransferFunction{T})
    return TransferFunction(ones(SISOTransferFunction{T}, t.outputs, t.inputs), t.inputs, t.outputs)
end

function one{T}(::Type{TransferFunction{T}})
    return TransferFunction(ones(SISOTransferFunction{T}, 1, 1), 1, 1)
end

## ZERO ##
zero(t::SISOTransferFunction) = SISOTransferFunction(zero(t.num), one(t.num))
zero{T}(::Type{SISOTransferFunction{T}}) = SISOTransferFunction(Poly([zero(T)], :s), Poly([one(T)], :s))
function zero{T}(t::TransferFunction{T})
    return TransferFunction(zeros(SISOTransferFunction{T}, t.outputs, t.inputs), t.inputs, t.outputs)
end

function zero{T}(::Type{TransferFunction{T}})
    return TransferFunction(zeros(SISOTransferFunction{T}, 1, 1), 1, 1)
end

## INDEXING ##
getindex(t::TransferFunction, i...) = t.matrix[i...] 
setindex!(t::TransferFunction, t2, i...) = (t.matrix[i...] = t2)
endof(t::TransferFunction) = length(t)

eltype{T}(::TransferFunction{T}) = T
eltype{T}(::SISOTransferFunction{T}) = T

function copy(t::TransferFunction)
    matrix = copy(t.matrix)
    return TransferFunction(matrix, t.inputs, t.outputs)
end

size(t::TransferFunction) = size(t.matrix)
length(t::TransferFunction) = length(t.matrix)

function gcd{T<:FloatingPoint, S<:FloatingPoint}(a::Poly{T}, b::Poly{S})
    #Finds the Greatest Common Denominator of two polynomials recursively using
    #Euclid's algorithm: http://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#Euclid.27s_algorithm
    if all(abs(b.a).<=2*eps(S))
        return a
    else
        s, r = divrem(a, b)
        return gcd(b, r)
    end
end

#####################################################################
##                         Math Operators                          ##
#####################################################################

## EQUALITY ##
==(t1::SISOTransferFunction, t2::SISOTransferFunction) = (t1.num == t2.num && t1.den == t2.den)

function ==(t1::TransferFunction, t2::TransferFunction)
    if size(t1) != size(t2)
        return false
    else
        return t1.matrix == t2.matrix
    end
end

## ADDITION ##
function +(t1::SISOTransferFunction, t2::SISOTransferFunction)
    num = t1.num*t2.den + t2.num*t1.den
    den = t1.den*t2.den
    return SISOTransferFunction(num, den)
end

function +(t::SISOTransferFunction, n::Real)
    num = t.num + n*t.den
    return SISOTransferFunction(num, t.den)
end
+(n::Real, t::SISOTransferFunction) = t + n

function +(t1::TransferFunction, t2::TransferFunction)
    #Check that the input-output sizes are consistent
    if t1.inputs != t2.inputs
        error("The first summand has %i input(s), but the second has %i.", t1.inputs, t2.inputs)
    end
    if t1.outputs != t2.outputs
        error("The first summand has %i output(s), but the second has %i.", t1.outputs, t2.outputs)
    end

    matrix = t1.matrix + t2.matrix
    return TransferFunction(matrix, t1.inputs, t1.outputs)
end

function +(t::TransferFunction, n::Real)
    if t.inputs == t.outputs == 1
        t2 = copy(t)
        t2.matrix[1,1] += n
    else
        error("Must be SISO to add a constant")
    end
    return t2
end
+(n::Real, t::TransferFunction) = t + n

## SUBTRACTION ##
function -(t1::SISOTransferFunction, t2::SISOTransferFunction)
    num = t1.num*t2.den - t2.num*t1.den
    den = t1.den*t2.den
    return SISOTransferFunction(num, den)
end

function -(n::Real, t::SISOTransferFunction)
    num = n*t.den - t.num
    return SISOTransferFunction(num, t.den)
end
-(t::SISOTransferFunction, n::Real) = t + -n

function -(t1::TransferFunction, t2::TransferFunction)
    #Check that the input-output sizes are consistent
    if t1.inputs != t2.inputs
        error("The first summand has %i input(s), but the second has %i.", t1.inputs, t2.inputs)
    end
    if t1.outputs != t2.outputs
        error("The first summand has %i output(s), but the second has %i.", t1.outputs, t2.outputs)
    end

    matrix = t1.matrix - t2.matrix
    return TransferFunction(matrix, t1.inputs, t1.outputs)
end

function -(n::Real, t::TransferFunction)
    if t.inputs == t.outputs == 1
        t2 = copy(t)
        t2.matrix[1,1] = n - t2.matrix[1,1]
    else
        error("Must be SISO to subtract a constant")
    end
    return t2
end
-(t::TransferFunction, n::Real) = t + -n

## NEGATION ##
-(t:: SISOTransferFunction) = SISOTransferFunction(-t.num, t.den)
-(t::TransferFunction) = TransferFunction(-t.matrix, t.inputs, t.outputs)

## MULTIPLICATION ##
function *(t1::SISOTransferFunction, t2::SISOTransferFunction)
    num = t1.num*t2.num
    den = t1.den*t2.den
    return SISOTransferFunction(num, den)
end

*(t::SISOTransferFunction, n::Real) = SISOTransferFunction(t.num*n, t.den)
*(n::Real, t::SISOTransferFunction) = t*n
.*(t::SISOTransferFunction, n::Real) = t*n
.*(n::Real, t::SISOTransferFunction) = t*n

function *(t1::TransferFunction, t2::TransferFunction)
    #Check that the input-output sizes are consistent
    if t1.inputs != t2.outputs
        error(@sprintf("Input->Output Mismatch: C = A*B: A has %i inputs, B has %i outputs", t1.inputs, t2.outputs))
    end

    matrix = t1.matrix * t2.matrix
    return TransferFunction(matrix, t1.inputs, t2.outputs)
end

function *(t::TransferFunction, n::Real)
    matrix = t.matrix*n
    return TransferFunction(matrix, t.inputs, t.outputs)
end
*(n::Real, t::TransferFunction) = t*n

## DIVISION ##
function /(t1::SISOTransferFunction, t2::SISOTransferFunction)
    num = t1.num*t2.den
    den = t1.den*t2.num
    return SISOTransferFunction(num, den)
end

function /(t1::TransferFunction, t2::TransferFunction)
    if t1.inputs == t2.inputs == t1.outputs == t2.outputs
        t = t1.matrix[1,1]/t2.matrix[1,1]
        matrix = reshape([t], 1, 1)
    else
        error("NotImplementedError: MIMO division isn't implemented")
    end
    return TransferFunction(matrix, t1.inputs, t2.outputs)
end

/{T<:Real}(t::TransferFunction, n::T) = t*(1/n)
/{T<:Real}(n::T, t::TransferFunction) = /(n*one(t), t)

#####################################################################
##                        Display Functions                        ##
#####################################################################
function print(io::IO, tf::SISOTransferFunction)
    #Convert the numerator and denominator to strings
    numstr = sprint(print, tf.num)
    denstr = sprint(print, tf.den)

    #Figure out the length of the separating line
    dashcount = max(length(numstr), length(denstr))

    #Center the numerator or denominator
    if length(numstr) < dashcount
        numstr = "$(repeat(" ", int((dashcount - length(numstr))/2)))$numstr"
    else
        denstr = "$(repeat(" ", int((dashcount - length(denstr))/2)))$denstr"
    end
    print(io, numstr)
    print(io, "\n$(repeat("-", dashcount))\n")
    print(io, denstr)
end

function print(io::IO, tf::TransferFunction)
    mimo = tf.inputs > 1 || tf.outputs > 1
    for i=1:tf.inputs
        for o=1:tf.outputs
            if mimo
                print(io, "Input $i to output $o\n")
            end
            print(io, tf.matrix[o, i])
        end
    end
end

function show(io::IO, tf::TransferFunction)
    print(io, "TransferFunction:\n")
    print(io, tf)
end
