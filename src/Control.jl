module Control
#Control Systems toolbox

export tf, TransferFunction

import Base: length, getindex, show, string, print
import Base: *, /, +, -


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
            jl_error("num and den output dimensions must match")
        end
        outputs = out_n

        if in_n != in_d
            jl_error("num and den input dimensions must match")
        end
        inputs = in_n

        ## Remove leading zeros on num and den ##
        #Data is deepcopied to prevent mutating the calling arrays
        data = {deepcopy(den), deepcopy(num)}
        for p in 1:2
            for o in 1:outputs
                for i in 1:inputs
                    nzfirst = 0
                    for j in 1:length(data[p][o,i])
                        #For each element in a num/den vector, check if the element is
                        #greater than ~0. This allows for small computation errors
                        #if input vector is calculated (i.e. 1e-17 ~0, and will be removed)
                        if abs(data[p][o,i][j]) > 2*eps(Float64)
                            nzfirst = j
                            break
                        end
                    end
                    if nzfirst == 0
                        # The array is all zeros
                        if p==1
                            #This is a denominator with zero value
                            jl_error("Input $i, output $o has a zero denominator")
                        else
                            #This is a numerator, make the denominator a 1
                            #Results in 0/1
                            data[2][o,i] = [0]
                            data[1][o,i] = [1]
                        end
                    else
                        #Truncate all leading zeros
                        data[p][o,i] = data[p][o,i][nzfirst:]
                    end
                end
            end
        end
        (new_den, new_num) = data
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

#TODO: Add simplification step, a + a - a != a currently, due to growing den
#TODO: Implement +, -, for scalars
#TODO: Implement *, /, neg

function +(self::TransferFunction, other::TransferFunction)
    ## Add two transfer function objects (parallel connection) ##

    #Check that the input-output sizes are consistent
    if self.inputs != other.inputs
        jl_error("The first summand has %i input(s), but the second has %i.", self.inputs, other.inputs)
    end
    if self.outputs != other.outputs
        jl_error("The first summand has %i output(s), but the second has %i.", self.outputs, other.outputs)
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

function -(self::TransferFunction, other::TransferFunction)
    ## Subtract two transfer function objects (parallel connection) ##

    #Check that the input-output sizes are consistent
    if self.inputs != other.inputs
        jl_error("The first summand has %i input(s), but the second has %i.", self.inputs, other.inputs)
    end
    if self.outputs != other.outputs
        jl_error("The first summand has %i output(s), but the second has %i.", self.outputs, other.outputs)
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

function _addSISO(num1::Vector{Float64}, den1::Vector{Float64}, num2::Vector{Float64}, den2::Vector{Float64})
    num = polyadd(polymul(num1, den2), polymul(num2, den1))
    den = polymul(den1, den2)
    return num, den
end

function polyadd(p1::Vector{Float64}, p2::Vector{Float64})
    n = length(p1)
    m = length(p2)
    if n > m
        p3 = Array(Float64, n)
        for i = 1:m
            p3[n-m+i] = p1[n-m+i] + p2[i]
        end
        for i = 1:n-m
            p3[i] = p1[i]
        end
    else
        p3 = Array(Float64, m)
        for i = 1:n
            p3[m-n+i] = p1[i] + p2[m-n+i]
        end
        for i = 1:m-n
            p3[i] = p2[i]
        end
    end
    return p3
end

function polymul(p1::Vector{Float64}, p2::Vector{Float64})
    n = length(p1)
    m = length(p2)
    if n == 0 || m == 0
        return Float64[]
    end
    p3 = zeros(Float64, n+m-1)
    for i = 1:length(p1)
        for j = 1:length(p2)
            p3[i+j-1] += p1[i] * p2[j]
        end
    end
    return p3
end

function polydiv(num::Vector{Float64}, den::Vector{Float64})
    m = length(den)
    n = length(num)

    if m == 0
        throw(DivideError())
    end

    deg = n-m+1
    if deg <= 0
        return Float64[0], num
    end
    d = zeros(Float64, n)
    q = zeros(Float64, deg)
    r = copy(num)
    for i = 1:deg
        quot = r[i] / den[1]
        q[i] = quot
        if i > 1
            d[i-1] = 0
            r[i-1] = 0
        end
        for j = 1:m
            k = i+j-1
            elem = den[j]*quot
            d[k] = elem
            r[k] -= elem
        end
    end
    r_mask = abs(r).>eps(Float64)
    if any(r_mask)
        return q, r[findfirst(r_mask):]
    else
        return q, r
    end
end

function gcd(a::Vector{Float64}, b::Vector{Float64})
    #Finds the Greatest Common Denominator of two polynomials using
    #Euclid's algorithm: http://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#Euclid.27s_algorithm
    
    if all(abs(b).<=2*eps(Float64))
        return a
    else
        s, r = polydiv(a, b)
        return gcd(b, r)
    end
end

function tf_simplify(num::Vector{Float64}, den::Vector{Float64})
    fact = gcd(num, den)
    if fact == [1]
        return (num, den)
    else
        num_s = polydiv(num, fact)[1]
        den_s = polydiv(den, fact)[1]
        return (num_s, den_s)
    end
end




#####################################################################
##                        Display Functions                        ##
#####################################################################

# TODO: Negative sign is dropped from first numerator
# TODO: Remove excess zeros on integers

function string(tf::TransferFunction)
    ## String representation of the transfer function ##

    mimo = tf.inputs > 1 || tf.outputs > 1
    outstr = ""

    for i=1:tf.inputs
        for j=1:tf.outputs
            if mimo
                outstr = "$outstrInput $i to output $j\n"
            end

            #Convert the numerator and denominator to strings
            numstr = _tfpolyToString(tf.num[j,i])
            denstr = _tfpolyToString(tf.den[j,i])

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

function _tfpolyToString{T}(p::Vector{T})
    #Convert polynomial vector into a string
    thestr = "0"

    #Compute the number of coefficients
    N = length(p)

    for k=1:(N)
        num = abs(p[k])
        if isa(num, Integer)
            coefstr = "$num"
        else
            coefstr = @sprintf("%.4f", num)
        end
        power = N-k
        if power == 0
            if coefstr != "0"
                newstr = coefstr
            else
                if k == 0
                    newstr = "0"
                else
                    newstr = ""
                end
            end
        elseif power == 1
            if coefstr == "0"
                newstr = ""
            elseif coefstr == "1"
                newstr = "s"
            else
                newstr = "$coefstr s"
            end
        else
            if coefstr == "0"
                newstr = ""
            elseif coefstr == "1"
                newstr = "s^$power"
            else
                newstr = "$coefstr s^$power"
            end
        end

        if k > 1
            if newstr != ""
                if p[k] < 0
                    thestr = "$thestr - $newstr"
                else
                    thestr = "$thestr + $newstr"
                end
            end
        elseif k == 0 && newstr != "" && p[k] < 0
            thestr = "-$newstr"
        else
            thestr = newstr
        end
    end
    return thestr
end

function show(io::IO, tf::TransferFunction)
    print("TransferFunction\n")
    print(string(tf))
end

#END MODULE
end
