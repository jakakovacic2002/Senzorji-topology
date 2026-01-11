using Mods, Primes, LinearAlgebra
using Random
using WGLMakie, ColorTypes

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

function maximal_complex_at_r(filtration, r_value)
    # Filter filtration to only include simplices with radius <= r_value
    complex_at_r = [simplex for (simplex, r) in filtration if r <= r_value]
    
    # Find maximal simplices
    max_simplices = find_maximal_simplices(complex_at_r)
    
    return max_simplices
end

function simplex_boundary(sx)
    # determine the codimension 1 simplices in the boundary of the simplex sx
    # 'remove' one entry of a tuple at each specific index by splicing the corresponding subtuples
    cd1faces = ([(sx[1:k-1]..., sx[k+1:end]...) for k in eachindex(sx)])
    return cd1faces
end

function simultaneous_reduce(A, B; findgens = false)
    # throw a 'wrong size' error
    if size(A)[2] != size(B)[1]
        throw(DimensionMismatch("cols(A) != rows(B)"))
    end
    # store the number of rows and columns of A 
    n_rows, n_cols = size(A)
    if findgens
        # initialize the basis of Cn (as the identity matrix of appropriate type)
        basis = Matrix{eltype(A)}(I, n_cols, n_cols)
    end
    # perform a column reduction on A, do the corresponding operations on columns of 'basis' and do the inverse operations on rows of B
    i, j = 1, 1
    last_nonzero_col = 0
    while i <= n_rows && j <= n_cols
        # check if A[i, j] can be used as a column pivot ...
        if A[i, j] == 0
            # if not, find a column pivot in the current row
            nonzero_col = j+1
            while nonzero_col <= n_cols && A[i, nonzero_col] == 0
                nonzero_col += 1;
            end
            # if there are only zeros in the current row to the right
            if nonzero_col == n_cols+1
                # go to a new row
                i += 1
                continue
            end
            # otherwise swap the columns of A ... 
            A[:, [j nonzero_col]] = A[:, [nonzero_col j]]
            # ... (and of 'basis') ...
            if findgens
                basis[:, [j nonzero_col]] = basis[:, [nonzero_col j]]
            end
            # ... and corresponding rows of B
            B[[j nonzero_col], :] = B[[nonzero_col j], :]
        end
        # store last nonzero column (for additional row reduction on B later)
        last_nonzero_col = j
        # ... and set A[i, j] as the new column pivot entry 
        pivot = A[i, j]
        # set that column pivot to 1 (in A) ...
        A[:, j] = A[:, j]/pivot
        # ... and do corresponding operation on 'basis' ...
        if findgens
            basis[:, j] = basis[:, j]/pivot
        end
        # and do the inverse operation on B
        B[j, :] = B[j, :]*pivot
        # annihilate the nonzero elements to the right of the pivot in A,
        # do the corresponding operations on 'basis'
        # and do the inverse operations on B
        for col in (j+1):n_cols
            scale = A[i, col]
            A[:, col] -= scale*A[:, j]
            if findgens
                basis[:, col] -= scale*basis[:, j] 
            end
            B[j, :] += scale*B[col, :]
        end
        # increase row and column indices and search for the next pivot
        i += 1
        j += 1
    end
    # now perform a row reduction on nonzero rows of B 
    # (Performing inverse operations on columns of A is not required. 
    # But we need to do the inverse operations on 'basis' to keep track of the generators.)
    n_rows, n_cols = size(B)
    i, j = last_nonzero_col+1, 1
    last_nonzero_row = last_nonzero_col
    while i <= n_rows && j <= n_cols
        # check if B[i, j] can be used as a column pivot ...
        if B[i, j] == 0
            # if not, find a row pivot in the current column
            nonzero_row = i+1
            while nonzero_row <= n_rows && B[nonzero_row, j] == 0
                nonzero_row += 1
            end
            # if there are only zeros in the current column downwards
            if nonzero_row == n_rows+1
                # go to a new column
                j += 1
                continue
            end
            # otherwise swap the rows of B ...
            B[[i nonzero_row], :] = B[[nonzero_row i], :]
            # ... and corresponding columns of 'basis'
            if findgens
                basis[:, [i nonzero_row]] = basis[:, [nonzero_row i]]
            end
            # (The corresponding columns of the column-reduced A are 0, no need to swap those...)
        end
        # store last nonzero row (to return the 'Betti number')
        last_nonzero_row = i
        # Set B[i, j] as the new pivot entry ...
        pivot = B[i, j]
        # ... and set that pivot to 1.
        B[i, :] = B[i, :]/pivot
        # Do the inverse operation on 'basis'.
        if findgens
            basis[:, i] = pivot*basis[:, i]
        end
        # annihilate the nonzero elements below the pivot in B
        # (No need for inverse operations on the columns of the reduced A...)
        for row in (i+1):n_rows
            scale = B[row, j]
            B[row, :] -= scale*B[i, :]
            # Do the inverse operation on 'basis'.
            if findgens
                basis[:, i] += scale*basis[:, row] 
            end
        end
        # increase row and column indices and search for the next pivot
        i += 1
        j += 1
    end
    if !findgens
        # return the corresponding Betti number ...
        return n_rows - last_nonzero_row
    else
        # ... or a named tuple (betti, gens) if findgens == true
        return (betti = n_rows - last_nonzero_row, gens = basis[:, (last_nonzero_row+1):end])
    end
end

function boundary_matrix(dimension_dict, n, p = 0)
    # check if p is 0 or prime
    if p != 0 && !isprime(p)
        throw(DomainError(p, "p should be 0 or a prime."))
    end
    # cases n = 0 and n = dim+1 should be treated separately
    if n == 0
        Dn = zeros(p==0 ? Float64 : Mod{p}, 1, length(dimension_dict[n]))
    elseif n == maximum(keys(dimension_dict))+1
        Dn = zeros(p==0 ? Float64 : Mod{p}, length(dimension_dict[n-1]), 1)
    else
        Dn = zeros(p==0 ? Float64 : Mod{p}, length(dimension_dict[n-1]), length(dimension_dict[n]))
        # otherwise, just replace the zeroes with 1's or -1's at appropriate places
        for j in eachindex(dimension_dict[n])
            boundary = simplex_boundary(dimension_dict[n][j])
            for i in eachindex(boundary)
                Dn[findfirst(==(boundary[i]), dimension_dict[n-1]), j] = -(-1)^i
            end
        end
    end
    return Dn
end

function dim_dict_scx(max_scx)
    # determine the dimension+1 of the complex
    n = maximum([length(sx) for sx in max_scx])
    # set up the dictionary (empty vectors of k-tuples in each dimension)
    dimension_dict = Dict((k-1) => Vector{NTuple{k, Int}}() for k in 1:n)
    # add all simplices in max_scx to the dictionary
    # (Sorting of vertices is necessary for some direct comparisons in a later function...)
    for sx in max_scx
        push!(dimension_dict[length(sx)-1], Tuple(sort(collect(sx))))
    end
    # fill the dictionary with all simplices from top down
    for k in (n-1):-1:1
        for sx in dimension_dict[k]
            union!(dimension_dict[k-1], [Tuple(sort(collect(x))) for x in simplex_boundary(sx)])
        end
    end
    return dimension_dict
end

function betti_numbers(max_scx, p; findgens = false)
    # build the 'dimension dictionary' of a simplicial complex given via its maximal simplices
    dimension_dict = dim_dict_scx(max_scx)
    # determine the dimension of max_scx
    d = maximum(keys(dimension_dict))
    # initialize the vector of Betti numbers ...
    betti = Int[]
    # ... and the vector of 'generators'
    if findgens
        gens = LinearCombination[]
    end
    # the current (first) boundary matrix
    Dk1 = boundary_matrix(dimension_dict, 0, p)
    # determine the Betti numbers (and generators) from two consecutive boundary matrices
    for k in 1:(d+1)
        # the next boundary matrix
        Dk = boundary_matrix(dimension_dict, k, p)
        # perform simultaneous reduction on D_{k-1} and D_k 
        # and store the corresponding Betti number (and 'generators')
        if findgens
            (β, genk) = simultaneous_reduce(Dk1, Dk; findgens = true)
            # store the Betti number
            push!(betti, β)
            # store each generator in dimension k-1 as a LinearCombination (of simplices of K)
            push!(gens, [LinearCombination(Dict(dimension_dict[k-1][scx_index] => genk[scx_index, col] for scx_index in eachindex(dimension_dict[k-1]))) for col in axes(genk, 2)]...)
        else
            push!(betti, simultaneous_reduce(Dk1, Dk; findgens = false))
        end
        # use Dk as the current (first) boundary matrix
        Dk1 = Dk
    end
    if !findgens
        # return a list of numbers if generators are not required
        return betti
    else
        # return a named tuple of two lists if generators are required
        return (betti = betti, gens = gens)
    end
end

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

function find_sphere_r(filtration, p=0)
    complex_at_r = Vector{Vector{Int}}()
    prev_r = -1.0
    
    for (simplex, r) in filtration
        # When we encounter a new r value, check the Betti numbers of the previous complex
        if r != prev_r && prev_r != -1.0 && !isempty(complex_at_r)
            # Find maximal simplices
            max_simplices = find_maximal_simplices(complex_at_r)
            
            # Convert to tuple format for betti_numbers
            max_scx = [Tuple(s) for s in max_simplices]
            
            # Calculate Betti numbers
            betti = betti_numbers(max_scx, p)
            
            # Check if Betti numbers are [1, 0, 1] (sphere topology)
            if length(betti) >= 3 && betti[1] == 1 && betti[2] == 0 && betti[3] == 1
                println("Found sphere topology (β₀=1, β₁=0, β₂=1) at R = $prev_r")
                return prev_r
            end
        end
        
        # Add current simplex to the complex
        push!(complex_at_r, simplex)
        prev_r = r
    end
    
    # Check the last r value
    if !isempty(complex_at_r)
        max_simplices = find_maximal_simplices(complex_at_r)
        max_scx = [Tuple(s) for s in max_simplices]
        betti = betti_numbers(max_scx, p)
        
        if length(betti) >= 3 && betti[1] == 1 && betti[2] == 0 && betti[3] == 1
            println("Found sphere topology (β₀=1, β₁=0, β₂=1) at R = $prev_r")
            return prev_r
        end
    end
    
    println("No complex with sphere topology (β₀=1, β₁=0, β₂=1) found")
    return nothing
end

# Helper function to find maximal simplices from a list of all simplices
function find_maximal_simplices(simplices)
    # Optimization: group simplices by dimension
    by_dim = Dict{Int, Vector{Vector{Int}}}()
    max_dim = 0
    
    for s in simplices
        dim = length(s)
        max_dim = max(max_dim, dim)
        if !haskey(by_dim, dim)
            by_dim[dim] = Vector{Vector{Int}}()
        end
        push!(by_dim[dim], s)
    end
    
    # Start with highest dimensional simplices (always maximal)
    maximal = Vector{Vector{Int}}()
    if haskey(by_dim, max_dim)
        append!(maximal, by_dim[max_dim])
    end
    
    # Check lower dimensions only if not contained in higher dimensions
    for dim in (max_dim-1):-1:1
        if !haskey(by_dim, dim)
            continue
        end
        
        for s1 in by_dim[dim]
            is_maximal = true
            # Only check against higher dimensional simplices
            for higher_dim in (dim+1):max_dim
                if !haskey(by_dim, higher_dim)
                    continue
                end
                for s2 in by_dim[higher_dim]
                    if all(v in s2 for v in s1)
                        is_maximal = false
                        break
                    end
                end
                if !is_maximal
                    break
                end
            end
            if is_maximal
                push!(maximal, s1)
            end
        end
    end
    
    return maximal
end

function plot_sphere_with_points_web(pts, edges, r)
    fig = Figure(size=(800, 800))
    ax = Axis3(fig[1,1], aspect=:data)

    # Sphere mesh
    θ = range(0, 2π, length=80)
    φ = range(0, π, length=40)
    x = [sin(phi)*cos(th) for phi in φ, th in θ]
    y = [sin(phi)*sin(th) for phi in φ, th in θ]
    z = [cos(phi) for phi in φ, th in θ]

    # Make sphere transparent
    # surface!(ax, x, y, z, color=RGBA(0.5, 0.7, 1.0, 0.2), transparency=true)

    # Extract coordinates from vector of vectors
    pts_x = [p[1] for p in pts]
    pts_y = [p[2] for p in pts]
    pts_z = [p[3] for p in pts]
    
    # Points
    scatter!(ax, pts_x, pts_y, pts_z, color=RGBA(1,0,0,1), markersize=15)

    # Draw filled circles around each point on the sphere
    for p in pts
        # Normalize to get point on unit sphere
        p_norm = p / norm(p)
        
        # Find two orthogonal vectors perpendicular to p_norm
        if abs(p_norm[3]) < 0.9
            v1 = cross(p_norm, [0.0, 0.0, 1.0])
        else
            v1 = cross(p_norm, [1.0, 0.0, 0.0])
        end
        v1 = v1 / norm(v1)
        v2 = cross(p_norm, v1)
        v2 = v2 / norm(v2)
        
        # Generate circle on sphere surface with radial segments for filling
        n_angles = 50
        angles = range(0, 2π, length=n_angles)
        n_radial = 20
        radial_steps = range(0, r, length=n_radial)
        
        # Create mesh points
        mesh_x = []
        mesh_y = []
        mesh_z = []
        
        for rad in radial_steps
            for angle in angles
                local_pt = cos(rad) * p_norm + sin(rad) * (cos(angle) * v1 + sin(angle) * v2)
                push!(mesh_x, local_pt[1])
                push!(mesh_y, local_pt[2])
                push!(mesh_z, local_pt[3])
            end
        end
        
        # Reshape for surface plot
        mesh_x_mat = reshape(mesh_x, n_angles, n_radial)
        mesh_y_mat = reshape(mesh_y, n_angles, n_radial)
        mesh_z_mat = reshape(mesh_z, n_angles, n_radial)
        
        surface!(ax, mesh_x_mat, mesh_y_mat, mesh_z_mat, color=RGBA(1.0, 0.5, 0.0, 1.0), shading=false)
    end

    # Edges
    for edge in edges
        i, j = edge[1], edge[2]
        lines!(ax, [pts[i][1], pts[j][1]], [pts[i][2], pts[j][2]], [pts[i][3], pts[j][3]], color=:black, linewidth=2)
    end

    fig
end

function plot_sphere_with_points(pts, r)
    fig = Figure(size=(800, 800))
    ax = Axis3(fig[1,1], aspect=:data)

    # Sphere mesh
    θ = range(0, 2π, length=80)
    φ = range(0, π, length=40)
    x = [sin(phi)*cos(th) for phi in φ, th in θ]
    y = [sin(phi)*sin(th) for phi in φ, th in θ]
    z = [cos(phi) for phi in φ, th in θ]

    # Make sphere transparent
    surface!(ax, x, y, z, color=RGBA(0.5, 0.7, 1.0, 0.2), transparency=true)

    # Extract coordinates from vector of vectors
    pts_x = [p[1] for p in pts]
    pts_y = [p[2] for p in pts]
    pts_z = [p[3] for p in pts]
    
    # Points
    scatter!(ax, pts_x, pts_y, pts_z, color=RGBA(1,0,0,1), markersize=15)

    # Draw filled circles around each point on the sphere
    for p in pts
        # Normalize to get point on unit sphere
        p_norm = p / norm(p)
        
        # Find two orthogonal vectors perpendicular to p_norm
        if abs(p_norm[3]) < 0.9
            v1 = cross(p_norm, [0.0, 0.0, 1.0])
        else
            v1 = cross(p_norm, [1.0, 0.0, 0.0])
        end
        v1 = v1 / norm(v1)
        v2 = cross(p_norm, v1)
        v2 = v2 / norm(v2)
        
        # Generate circle on sphere surface with radial segments for filling
        n_angles = 50
        angles = range(0, 2π, length=n_angles)
        n_radial = 20
        radial_steps = range(0, r, length=n_radial)
        
        # Create mesh points
        mesh_x = []
        mesh_y = []
        mesh_z = []
        
        for rad in radial_steps
            for angle in angles
                local_pt = cos(rad) * p_norm + sin(rad) * (cos(angle) * v1 + sin(angle) * v2)
                push!(mesh_x, local_pt[1])
                push!(mesh_y, local_pt[2])
                push!(mesh_z, local_pt[3])
            end
        end
        
        # Reshape for surface plot
        mesh_x_mat = reshape(mesh_x, n_angles, n_radial)
        mesh_y_mat = reshape(mesh_y, n_angles, n_radial)
        mesh_z_mat = reshape(mesh_z, n_angles, n_radial)
        
        surface!(ax, mesh_x_mat, mesh_y_mat, mesh_z_mat, color=RGBA(1.0, 0.5, 0.0, 1.0), transparency=true, shading=false)
        
        # Draw black outline at the boundary of the circle
        outline_x = Float64[]
        outline_y = Float64[]
        outline_z = Float64[]
        for angle in angles
            local_pt = cos(r) * p_norm + sin(r) * (cos(angle) * v1 + sin(angle) * v2)
            push!(outline_x, local_pt[1])
            push!(outline_y, local_pt[2])
            push!(outline_z, local_pt[3])
        end
        lines!(ax, outline_x, outline_y, outline_z, color=:black, linewidth=2)
    end

    fig
end

# REMOVE POINTS

function find_max_redundant_points_maximal(max_complex, n_points, p=0)
    println("Starting with $(length(max_complex)) maximal simplices")
    
    redundant_sets = Vector{Vector{Int}}[]
    
    # Try removing k points for k = 1, 2, 3, ...
    for k in 1:(n_points-4)  # Need at least 4 points for a 2-sphere
        println("Trying to remove $k points...")
        
        # Generate all combinations of k points to remove
        vertices = collect(1:n_points)
        
        # Collect all valid combinations for this k
        valid_combinations = collect_valid_combinations(max_complex, vertices, k, p)
        
        if !isempty(valid_combinations)
            redundant_sets = valid_combinations
            println("✓ Can remove $k points while maintaining sphere topology")
            println("  Found $(length(valid_combinations)) valid combination(s)")
        else
            println("✗ Cannot remove $k points")
            break  # If we can't remove k points, we can't remove k+1 either
        end
    end
    
    return redundant_sets
end

function collect_valid_combinations(max_complex, vertices, k, p)
    valid_combos = Vector{Vector{Int}}()
    collect_combinations_recursive(max_complex, vertices, k, Int[], 1, p, valid_combos)
    return valid_combos
end

function collect_combinations_recursive(max_complex, vertices, k, removed, start_idx, p, valid_combos)
    # Check topology after each removal
    if !isempty(removed)
        # Build new maximal complex directly
        new_max_complex = Vector{Vector{Int}}()
        boundary_faces = Vector{Vector{Int}}()
        
        # Separate simplices: kept vs deleted
        for simplex in max_complex
            if any(v in removed for v in simplex)
                # This simplex contains a removed vertex - generate boundary faces
                for i in 1:length(simplex)
                    face = [simplex[1:i-1]..., simplex[i+1:end]...]
                    # Only consider faces that don't contain removed vertices
                    if !isempty(face) && !any(v in removed for v in face)
                        if !(face in boundary_faces)
                            push!(boundary_faces, face)
                        end
                    end
                end
            else
                # Keep this simplex (it's still maximal)
                push!(new_max_complex, simplex)
            end
        end
        
        # Add boundary faces that are NOT subsets of kept simplices
        for face in boundary_faces
            is_maximal = true
            for kept_simplex in new_max_complex
                # Check if face is a proper subset of kept_simplex
                if length(face) < length(kept_simplex) && all(v in kept_simplex for v in face)
                    is_maximal = false
                    break
                end
            end
            if is_maximal
                push!(new_max_complex, face)
            end
        end
        
        if isempty(new_max_complex)
            return  # Topology broken, stop this branch
        end
        
        # Convert to tuple format and check Betti numbers
        max_scx = [Tuple(s) for s in new_max_complex]
        betti = betti_numbers(max_scx, p)
        
        # Check if it's still a sphere (β₀=1, β₁=0, β₂=1)
        if !(length(betti) >= 3 && betti[1] == 1 && betti[2] == 0 && betti[3] == 1)
            return  # Topology broken, stop this branch
        end
        
        # If we've removed k points and topology is still good, save this combination!
        if length(removed) == k
            push!(valid_combos, copy(removed))
            return
        end
    end
    
    # Base case: we've successfully removed k points
    if length(removed) == k
        return
    end
    
    # Recursive case: try adding each remaining vertex to removed set
    for i in start_idx:length(vertices)
        push!(removed, vertices[i])
        collect_combinations_recursive(max_complex, vertices, k, removed, i+1, p, valid_combos)
        pop!(removed)
    end
end

# Call this to get the minimal R value
function get_r(points)
    filtration = cech_filtration(points)
    return find_sphere_r(filtration)
end
