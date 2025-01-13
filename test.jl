struct Test
    a::Int64
    b::Int64
end

module inner
using Main: Test
function func(test::Test)
    println(plus(test.a, test.b))    
end

function plus(a::Int64, b::Int64)   
    return a + b
end

export func
end

test = Test(1, 2)
func(test)

using .inner


