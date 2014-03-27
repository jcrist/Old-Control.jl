using Base.Test
include("../src/Control.jl")
using Control

#############################################################################
# Set up Custom Test Framework
type TestData
    npass::Integer
    nfail::Integer
    nerr::Integer
    fail::Vector{Test.Failure}
    error::Vector{Test.Error}
end
TestData() = TestData(0, 0, 0, Test.Failure[], Test.Error[])

function custom_handler(r::Test.Success)
    print("\033[92m.\033[0m")
    TEST_DATA.npass += 1
end

function custom_handler(r::Test.Failure)
    print("\033[91mF\033[0m")
    TEST_DATA.nfail += 1
    append!(TEST_DATA.fail, [r])
end

function custom_handler(r::Test.Error)
    print("\033[93mE\033[0m")
    TEST_DATA.nerr += 1
    append!(TEST_DATA.error, [r])
end

info_string = "Julia Version: $(Base.VERSION), $(Sys.MACHINE)"

##############################################################################

# List of tests
path = Base.splitdir(Base.source_path())[1]
tests = [splitext(x)[1] for x in filter(x->contains(x, "test_"), readdir(path))]

# Initialize global test data struct
global TEST_DATA = TestData()

println(info_string)
# Run the tests
println("Running Tests:")
Test.with_handler(custom_handler) do
    for t in tests
        print(" * $(t): ")
        include("$(t).jl")
        print("\n")
    end
end
print("\n")

# Print failure and error information
if TEST_DATA.nfail + TEST_DATA.nerr > 0
    println("$(repeat("-", 80))")
    for r in TEST_DATA.fail
        println("\033[91mTest failed in:\033[0m $(r.expr)\n")
    end
    for r in TEST_DATA.error
        println("\033[93mTest errored in:\033[0m $(r.expr)")
        println(r.err)
        Base.show_backtrace(STDOUT, :do_test, r.backtrace, 1:typemax(Int))
        println("\n")
    end
end

# Print test summary
println(repeat("-", 80))
println("Passed: $(TEST_DATA.npass), Failed: $(TEST_DATA.nfail), Errored: $(TEST_DATA.nerr)")
println(repeat("-", 80))

# Print if ok, in case above wasn't clear enough ;)
if TEST_DATA.nfail + TEST_DATA.nerr == 0
    println("OK TO COMMIT\n")
else
    println("DO NOT COMMIT\n")
end
