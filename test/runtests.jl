using Base.Test
include("../src/Control.jl")
using Control

tests = ["transferfunction",
         "statespace"]

println("Running Tests:")
for t in tests
    println(" * $(t)")
    include("$(t).jl")
end
