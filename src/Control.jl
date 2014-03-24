#Module for working with control systems

#Designed to work similiar to Matlab or Python Control Systems toolbox

#Submodules:
# *transferfunction:
#      Transfer Function type definition, and operator overloading
# *statespace:
#      State Space type definition, and operator overloading

module Control

export TransferFunction, tf, ss2tf
export StateSpace, ss

import Base: length, getindex, show, string, print, size
import Base: *, /, +, -, ==

#Abstract type for all system types
abstract Sys

include("polylib.jl")
include("statespace.jl")
include("transferfunction.jl")


end     #Module
