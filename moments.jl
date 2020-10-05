using Base
using LinearAlgebra
using FileIO
using GeometryBasics


@doc """ because julia defaults are dumb"""->
function dropsum(tensor,dims)
    dropdims(sum(tensor,dims=dims),dims=(isa(dims,Integer) ? (dims,) : dims))
end

@doc """ computes the volume of an obj from vertices of the boundary mesh"""->
function vol(mesh_verts)
    det(mesh_verts)/6
end

@doc """ (n,d) -> (d)"""->
function com(mesh_verts)
    dropsum(mesh_verts,1)/4
end

@doc """ (n,d), (d) -> (d,d)"""->
function ExxT(V,mu)
    (V'*V)/20+(4/5)*mu*mu'
end

@doc """ Given the coordinates for the triangle vertices of a mesh (N,n,d)
        computes the volume, center of mass (d), and covariance matrix (d,d)"""->
function compute_moments(triangles)
    ws = [vol(tri) for tri = triangles]
    Vol = sum(ws)
    ws /= Vol
    coms = vcat([com(tri)' for tri = triangles]...) # (N,3)
    mu = dropsum(coms.*ws,1)
    xxT = sum([ExxT(tri,com)*w for (tri,com,w)=zip(triangles,eachrow(coms),ws)])
    covar = xxT-mu*mu'
    Vol,mu,covar
end

@doc """ converts obtuse GeometryBasics.jl abstract triangles to 3x3 arrays"""->
function notDumbTriangle(triangle)
    Array{Float32,2}(hcat(triangle.points.data...))'
end

# Test it out:
mesh = load("cube.obj")
tris = [notDumbTriangle(tri) for tri=mesh]
v1,m1,C1 = compute_moments(tris)
offset =[1 -2 .5]
v2,m2,C2 = compute_moments([tr.+offset for tr=tris])
volerr = abs(v1-v2)/abs(v1+v2)
meanerr = norm(m2-offset' - m1)/norm(m2+m1)
covarerr = norm(C1-C2)/norm(C1+C2)
println("Relative errors in V:$volerr,Com:$meanerr,Covar:$covarerr")