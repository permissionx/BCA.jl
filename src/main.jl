using LinearAlgebra
import Base: push!
include("types.jl")
include("bca.jl")  # In namespace BCA: BCA->(QLoss, Constants)
include("io.jl")
include("geometry.jl")
include("dynamics.jl")


# modules


#simulator, atom_p, atom_t = test_collision()


# 近似：
# 1. 多个碰撞一个
# 2. 间隙子位移能
# 3. 替代位位移能

# 未解决：
# 1. 入射粒子的位移能
# 2. 间隙子的坐标

