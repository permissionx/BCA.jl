mutable struct Son
    p::Vector{Int64}
end

mutable struct Father
    sons::Vector{Son}
end




p = [1,2,3]
son = Son(p)
p = [4,5,6]
son2 = Son(p2)
father_ori = Father([son,son2])
father = deepcopy(father_ori)
father.sons = copy(father_ori.sons)
father.sons[1].p = [100,200,300]
println(father_ori.sons[1].p)


