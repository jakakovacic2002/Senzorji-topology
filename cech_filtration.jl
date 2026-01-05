using LinearAlgebra
using Random
using LightGraphs

struct Ball
    center::Vector{Float64}
    radius::Float64
end

function circumcircle(p1::Vector{Float64},
                      p2::Vector{Float64},
                      p3::Vector{Float64})

    a = p2 - p1
    b = p3 - p1

    # Degeneracy check
    if norm(cross(a, b)) < 1e-12
        error("Points are collinear; circumcircle undefined")
    end

    G = [
        dot(a, a)  dot(a, b)
        dot(a, b)  dot(b, b)
    ]

    rhs = 0.5 * [
        dot(a, a)
        dot(b, b)
    ]

    coeffs = G \ rhs
    center = p1 + coeffs[1] * a + coeffs[2] * b
    radius = norm(center - p1)

    return Ball(center, radius)
end

function circumsphere(p1::Vector{Float64},
                      p2::Vector{Float64},
                      p3::Vector{Float64},
                      p4::Vector{Float64})

    a = p2 - p1
    b = p3 - p1
    c = p4 - p1

    # Check degeneracy
    if abs(det(hcat(a, b, c))) < 1e-12
        error("Points are coplanar; circumsphere undefined")
    end

    A = 2 * hcat(a, b, c)'
    rhs = [
        dot(p2,p2) - dot(p1,p1),
        dot(p3,p3) - dot(p1,p1),
        dot(p4,p4) - dot(p1,p1)
    ]

    center = A \ rhs
    radius = norm(center - p1)

    return Ball(center, radius)
end

function ball_from_boundary(n::Vector{Vector{Float64}})
    k = length(n)

    if k == 0
        return Ball(zeros(3), 0.0)

    elseif k == 1
        return Ball(n[1], 0.0)

    elseif k == 2
        c = (n[1] + n[2]) / 2
        return Ball(c, norm(n[1] - c))

    elseif k == 3
        return circumcircle(n[1], n[2], n[3])

    elseif k == 4
        return circumsphere(n[1], n[2], n[3], n[4])

    else
        error("Boundary has more than d+1 points")
    end
end

function miniball(t::Vector{Vector{Float64}},
                  n::Vector{Vector{Float64}})
    # Base case: no more points to process
    if isempty(t) || length(n) == 4
        return ball_from_boundary(n)
    end

    # Choose u ∈ t
    u = pop!(t)

    # Recursive call without u
    B = miniball(t, n)

    # If u is already inside the ball, nothing changes
    if norm(u - B.center) ≤ B.radius + 1e-10
        push!(t, u)
        return B
    end

    # Otherwise, u must lie on the boundary
    push!(n, u)
    B2 = miniball(t, n)
    pop!(n)
    push!(t, u)

    return B2
end

function minimum_enclosing_ball(s::Vector{Vector{Float64}})
    t = copy(s)
    shuffle!(t)
    n = Vector{Vector{Float64}}()
    return miniball(t, n)
end

function pairwise_distances(points)
    n = length(points)
    D = zeros(n, n)
    for i in 1:n, j in i+1:n
        d = norm(points[i] - points[j])
        D[i,j] = d
        D[j,i] = d
    end
    return D
end

function adjacency_graph(points, r)
    n = length(points)
    adj = [BitSet() for _ in 1:n]

    for i in 1:n
        for j in i+1:n
            if norm(points[i] - points[j]) <= 2r
                push!(adj[i], j)
                push!(adj[j], i)
            end
        end
    end
    return adj
end

function cech_filtration(points::Vector{Vector{Float64}})
    D = pairwise_distances(points)
    filtration = []
    n = length(points)

    # Vertices
    for i in 1:n
        push!(filtration, ([i], 0.0))
    end

    # Edges
    for i in 1:n
        for j in i+1:n
            r = D[i,j] / 2
            push!(filtration, ([i, j], r))
        end
    end

    # Triangles
    for i in 1:n
        for j in i+1:n
            for k in j+1:n
                pts = [points[i], points[j], points[k]]
                B = minimum_enclosing_ball(pts)
                push!(filtration, ([i, j, k], B.radius))
            end
        end
    end

    # Tetrahedra
    for i in 1:n
        for j in i+1:n
            for k in j+1:n
                for l in k+1:n
                    pts = [points[i], points[j], points[k], points[l]]
                    B = minimum_enclosing_ball(pts)
                    push!(filtration, ([i, j, k, l], B.radius))
                end
            end
        end
    end

    sort!(filtration, by = x -> x[2])

    return filtration
end

# function r_complex(filtration, r)
#     complex = Vector{Vector{Int}}()
#     for (simplex, r_s) in filtration
#         if r_s <= r
#             push!(complex, simplex)
#         end
#     end

#     return complex
# end

# function is_connected_complex(complex::Vector{Vector{Int}})
#     # collect vertices
#     vertices = Set{Int}()
#     for simplex in complex
#         for v in simplex
#             push!(vertices, v)
#         end
#     end

#     if isempty(vertices)
#         return false
#     end

#     # adjacency list
#     adj = Dict{Int, Set{Int}}()
#     for v in vertices
#         adj[v] = Set{Int}()
#     end

#     # add edges (1-simplices)
#     for simplex in complex
#         if length(simplex) == 2
#             a, b = simplex
#             push!(adj[a], b)
#             push!(adj[b], a)
#         end
#     end

#     # DFS
#     visited = Set{Int}()
#     stack = [first(vertices)]

#     while !isempty(stack)
#         v = pop!(stack)
#         if v ∉ visited
#             push!(visited, v)
#             for u in adj[v]
#                 if u ∉ visited
#                     push!(stack, u)
#                 end
#             end
#         end
#     end

#     return length(visited) == length(vertices)
# end

function smallest_r(filtration)
    χ = 0
    prev_r = -1.0
    r_values = Float64[]
    
    for (simplex, r) in filtration
        if r != prev_r
            if prev_r != -1.0 && χ == 2
                println("Found Euler characteristic = 2 at R = $prev_r")
                push!(r_values, prev_r)
            end
            prev_r = r
        end
        
        k = length(simplex) - 1
        χ += (-1)^k
    end
    
    if χ == 2
        println("Found Euler characteristic = 2 at R = $prev_r")
        push!(r_values, prev_r)
    end
    
    if isempty(r_values)
        println("No complex with Euler characteristic = 2 found")
    end
    
    return r_values
end

