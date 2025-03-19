mutable struct A
    a::Vector{Int64}
end

mutable struct B
    a::Vector{Int64}
end

a = A([1,2,3])
b = B(a.a)

