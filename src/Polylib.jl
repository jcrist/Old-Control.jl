module Polylib
#Module for dealing with polynomials

#NOTE: Most of this code is a modified version of the Polynomial.jl library
# found at https://github.com/vtjnash/Polynomial.jl. I didn't like how this was
# done, using a special type as a polynomial, instead of just an 1d array. This
# is a forked version for use in this REPO, doing polynomials similiar to how
# Numpy or Matlab does them. Some functions have also been added, in addition to
# the functions found in vtjnash's version.

#Currently, all polynomials must be of type Vector{Float64}. This is because
#the control.jl module uses this type extensively. In the future, these could
#easily be paramatrized to form a more general library.

#TODO:
#--> Paramatrize all functions to allow for different types
#--> Documentation
#--> Unit Tests

export polyadd, polymul, polydiv, polygcd, polyfracsimp, polytostring

import Base: length, getindex, show, string, print

function polyadd(p1::Vector{Float64}, p2::Vector{Float64})
    #Adds two polynomials.
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
    #Multiplies two polynomials
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
    #Divides two polynomials. Returns (Solution, Remainder).
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

function polygcd(a::Vector{Float64}, b::Vector{Float64})
    #Finds the Greatest Common Denominator of two polynomials recursively using
    #Euclid's algorithm: http://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#Euclid.27s_algorithm
    #TODO: Rewrite in procedural form
    
    if all(abs(b).<=2*eps(Float64))
        return a
    else
        s, r = polydiv(a, b)
        return polygcd(b, r)
    end
end

function polyfracsimp(num::Vector{Float64}, den::Vector{Float64})
    #Simpifies a polynomial fraction by factoring out the greatest common
    #denominator from both the numerator and denominator
    fact = polygcd(num, den)
    if fact == [1]
        return (num, den)
    else
        num_s = polydiv(num, fact)[1]
        den_s = polydiv(den, fact)[1]
        return (num_s, den_s)
    end
end

function polytostring{T}(p::Vector{T}; var="s")
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
                newstr = var
            else
                newstr = "$coefstr $var"
            end
        else
            if coefstr == "0"
                newstr = ""
            elseif coefstr == "1"
                newstr = "$var^$power"
            else
                newstr = "$coefstr $var^$power"
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

end     #Module
