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
    print_with_color(:green, ".")
    TEST_DATA.npass += 1
end

function custom_handler(r::Test.Failure)
    print_with_color(:red, "F")
    TEST_DATA.nfail += 1
    append!(TEST_DATA.fail, [r])
end

function custom_handler(r::Test.Error)
    print_with_color(:yellow, "E")
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

println(info_string, "\n")
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
        print_with_color(:red, "Test failed in: ")
        println("$(r.expr)\n")
    end
    for r in TEST_DATA.error
        print_with_color(:yellow, "Test errored in: ")
        println("$(r.expr)")
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
    print_with_color(:green, "[OK]\n\n")
    exit(0)
else
    print_with_color(:red, "[Fail]\n\n")
    exit(1)
end
