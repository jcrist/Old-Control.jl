#Module for working with control systems

#Designed to work similiar to Matlab or Python Control Systems toolbox

#Submodules:
# *transferfunction:
#      Transfer Function type definition, and operator overloading
# *statespace:
#      State Space type definition, and operator overloading

module Control

export TransferFunction, tf
export StateSpace, ss

import Base: length, size, one, zero, copy, convert, promote_rule
import Base: getindex, setindex!, eltype, show, print 
import Base: *, /, +, -, ==, .*
using Polynomial

#Try to import slicot
if haskey(Pkg.installed(), "Slicot")
    import Slicot.Simple: tb04ad
    SLICOT_DEFINED = true
else
    SLICOT_DEFINED = false
end

if SLICOT_DEFINED
    export ss2tf
end

#Abstract type for all system types
abstract Sys

#Major type definitions
include("statespace.jl")
include("transferfunction.jl")

end     #Module
