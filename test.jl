#using Revise
includet("constraints.jl")
using .constraints
#constraint
#r = rand_scene(2)

using ForwardDiff

function combined(x)
    constraint(state2scene(x))
end
#combined(rand(4,3,1))
g = x -> ForwardDiff.jacobian(combined,x)

g(rand(4,3,1))

# g(rand(4,3,3))