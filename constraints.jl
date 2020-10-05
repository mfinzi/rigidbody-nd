
# using ForwardDiff

# f(x::Vector) = x;

# x = rand(5)

# g = x -> ForwardDiff.jacobian(f, x); # g = Df

# g(x)

# body
# / \ \
#rigid articulated soft

# Rigid Body
#   - mesh
#   - n-pos

# Articulated Body
#  [Rigid body parts]

# scene
# list of [Bodies] for now

# func: state(scene) -> x,v as Vector
module constraints
export Body,RigidBody,Scene,constraint,rand_scene,state2scene

using Base
using LinearAlgebra

abstract type Body end

mutable struct RigidBody{T} <: Body
    X::Array{T,2} #(4,3)
    moments::Array{Float64,1}
    #mesh::
end

struct Scene
    bodies::Array{Body,1}
end


function constraint(s::Scene)
    vcat([constraint(body) for body=s.bodies]...)
end

function constraint(body::RigidBody{T}) where T
    Phi = Array{T,1}(undef,6)
    k =1
    for i in 1:4
        for j in 1:i-1
            Phi[k] = sum((body.X[i,:]-body.X[j,:]).^2)
            k+=1
        end
    end
    Phi
end


function rand_scene(N::Integer)
    Scene([RigidBody(rand(4,3),rand(4)) for i=1:N])
end

function state(s::Scene)
    cat(3,[state(body) for body=s.bodies]...)
end

function state(body::RigidBody)
    body.X
end

function state2scene(X)
    N = size(X)[3]
    Scene([RigidBody(X[:,:,i],rand(4)) for i=1:N])
end


function \delta 

# function inject_state!(s::Scene,X)
#     for (i,body) in enumerate(s.bodies)
#         body.X = X[:,:,i]
#     end
# end

end