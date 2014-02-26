#Module for working with control systems

#Designed to work similiar to Matlab or Python Control Systems toolbox

#Submodules:
# *transferfunction:
#      Transfer Function type definition, and operator overloading
# *statespace:
#      State Space type definition, and operator overloading
# *polylib:
#      Alternative to Polynomial.jl, accomplishes what I need it to do

module Control

#Abstract type for all system types
abstract Sys

import transferfunction: TransferFunction, tf, ss2tf
export TransferFunction, tf, ss2tf

import statespace: StateSpace, ss
export StateSpace, ss

end     #Module
