module Control

export tf, TransferFunction

import Base: length, getindex, show, string, print
import Base: *, /, +, -
using Polynomial 

type TransferFunction
    num::Poly
    den::Poly
end

function tf(num::Array, den::Array)
    TransferFunction(Poly(num), Poly(den))
end

*(a::TransferFunction, b::TransferFunction) = TransferFunction(a.num*b.num, a.den*b.den)

function poly_to_string(p::Poly)
    #Convert Polynomial to string
    n = length(p)
    if n <= 0
        result = "0"
    elseif n == 1
        result = "$(p[1])"
    else
        result = "$(p[1])s^$(n-1)"
        for i = 2:n-1
            if p[i] != 0
                result = "$result + $(p[i])s^$(n-i)"
            end
        end
        if p[n] != 0
            result = "$result + $(p[n])"
        end
    end
    return result
end 

function string(tf::TransferFunction)
    n_str = poly_to_string(tf.num)
    d_str = poly_to_string(tf.den)
    n = max(length(n_str), length(d_str))
    r = _center(n_str, n)
    r = "$r\n$(repeat("-", n))\n"
    r = "$r$(_center(d_str, n))"
    return r
end

function show(io::IO, tf::TransferFunction)
    print(string(tf))
end
    
function _center(str::String, w::Int)
    lbound = ifloor((w - length(str))/2)
    buffer = repeat(" ", lbound)
    "$(buffer)$(str)$(buffer)"
end
    
end
