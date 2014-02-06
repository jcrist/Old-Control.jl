#Module for working with control systems

#Designed to work similiar to Matlab or Python Control Systems toolbox

#Submodules:
#-->XferFcn:  Transfer Function type definition, and operator overloading
#-->Polylib:  Alternative to Polynomial.jl, accomplishes what I need it to do

module Control

import XferFcn: TransferFunction, tf
export TransferFunction, tf

end     #Module
