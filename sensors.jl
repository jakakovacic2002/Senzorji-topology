using DelimitedFiles, LinearAlgebra, WGLMakie, ColorTypes, BoundingSphere, Mods, Primes, Random

function read_points(filename)
    raw = read(filename, String)

    # Remove the enclosing double braces {{ ... }}
    inner = raw[3:end-2]          # drop first two chars "{{" and last "}}"

    # Split into individual vector strings
    vecs = split(inner, "},{")

    # Parse each vector
    points = [parse.(Float64, split(replace(v, ['{', '}'] => ""), ",")) for v in vecs]

    # Convert to a matrix (Nx3)
    return reduce(vcat, (p' for p in points))
end

#VR CODE
function distance_table(pts)
    n = size(pts, 1)
    D = zeros(Float64, n, n)
    for i in 1:(n-1)
        for j in (i+1):n
            D[i, j] = D[j, i] = norm(pts[i, :] - pts[j, :])
        end
    end
    return D
end
function VR(r, D)
    edges = Tuple{Int,Int}[]
    n = size(D, 1)
    for i in 1:(n-1)
        for j in (i+1):n
            if D[i, j] ≤ r
                push!(edges, (i, j))
            end
        end
    end
    return edges
end
#my hw1 code for connected components
function edgedictionary(V, E)
    n = V[end]
    # set up the dictionary (empty sets of k-tuples in each dimension)
    edge_dict = Dict(k => Set{Int}() for k in 1:n)
    for (v1, v2) in E
        push!(edge_dict[v1], v2)
        push!(edge_dict[v2], v1)
    end
    return edge_dict
end
function component_of(edge_dict, vertex, visited=Set{Int}())
    push!(visited, vertex)
    component = [vertex]

    for neighbor in edge_dict[vertex]
        if !(neighbor in visited)
            append!(component, component_of(edge_dict, neighbor, visited))
        end
    end
    return sort(component)
end
function findcomponents(V, E)
    edge_dict = edgedictionary(V, E)

    visited = Set{Int}()
    components = []

    for vertex in V
        if !(vertex in visited)
            push!(components, component_of(edge_dict, vertex, visited))
        end
    end
    return components
end

function is_connected(V, E; exclude=Set{Int}())
    if isempty(exclude)
        comp = findcomponents(V, E)
    else
        E = [e for e in E if !(e[1] in exclude || e[2] in exclude)]
        comp = findcomponents(V, E)
        comp = [c for c in comp if all(x -> !(x in exclude), c)]
    end
    return (length(comp) == 1)
end
function find_r(pts; steps=1000)
    D = distance_table(pts)
    V = collect(1:(size(pts, 1)))
    lo = 0.0
    hi = maximum(D)
    for _ in 1:steps
        r = (lo + hi) / 2
        edges = VR(r, D)
        if is_connected(V, edges)
            hi = r
        else
            lo = r
        end
    end
    return (lo + hi) / 2
end

function plot_edges(pts, edges)
    fig = Figure(size=(800, 800))
    ax = Axis3(fig[1, 1], aspect=:data)

    # Sphere mesh
    θ = range(0, 2π, length=80)
    φ = range(0, π, length=40)
    x = [sin(phi) * cos(th) for phi in φ, th in θ]
    y = [sin(phi) * sin(th) for phi in φ, th in θ]
    z = [cos(phi) for phi in φ, th in θ]

    # Make sphere transparent
    surface!(ax, x, y, z, color=RGBA(0.5, 0.7, 1.0, 0.2), transparency=true)

    # Points
    scatter!(ax, pts[:, 1], pts[:, 2], pts[:, 3], color=RGBA(1, 0, 0, 1), markersize=15)

    # Edges
    for (i, j) in edges
        lines!(ax, [pts[i, 1], pts[j, 1]], [pts[i, 2], pts[j, 2]], [pts[i, 3], pts[j, 3]], color=:black, linewidth=2)
    end

    display(fig)
end




##############
#
# This is for VR conectivity on DATASET
#
##############
function r_preset()
    pts = read_points("A_sensors_data.txt")
    r = find_r(pts)
    D = distance_table(pts)
    #V = collect(1:(size(pts,1)))
    E = VR(r, D)
    plot_edges(pts, E)

    return r
end




# ČECH CODE
function MiniBall_edges(points)
    n = size(points, 1)
    edges = Dict{Tuple{Int,Int},Float64}()
    for i in 1:n-1, j in i+1:n
        edges[(i, j)] = norm(points[i, :] - points[j, :]) / 2
    end
    return edges
end
function MiniBall_triangles(points)
    n = size(points, 1)
    tris = Dict{NTuple{3,Int},Float64}()

    for i in 1:n-2, j in i+1:n-1, k in j+1:n
        tris[(i, j, k)] = boundingsphere([
            Vector(points[i, :]),
            Vector(points[j, :]),
            Vector(points[k, :])
        ])[2]
    end
    return tris
end
function MiniBall_tetrahedra(points)
    n = size(points, 1)
    tets = Dict{NTuple{4,Int},Float64}()

    for i in 1:n-3, j in i+1:n-2, k in j+1:n-1, l in k+1:n
        tets[(i, j, k, l)] = boundingsphere([
            Vector(points[i, :]),
            Vector(points[j, :]),
            Vector(points[k, :]),
            Vector(points[l, :])
        ])[2]
    end
    return tets
end
#cech already returns a dim_dict
function cech(R, points, MB_edges, MB_triangles, MB_tetrahedra)
    n = size(points, 1)
    simplices = Dict((k - 1) => Vector{NTuple{k,Int}}() for k in 1:4)
    simplices[0] = [(i,) for i in 1:n]
    simplices[1] = [(i, j) for ((i, j), r) in MB_edges if r ≤ R]
    simplices[2] = [(i, j, k) for ((i, j, k), r) in MB_triangles if r ≤ R]
    simplices[3] = [(i, j, k, l) for ((i, j, k, l), r) in MB_tetrahedra if r ≤ R]
    return simplices
end
function simplex_boundary(sx)
    # determine the codimension 1 simplices in the boundary of the simplex sx
    # 'remove' one entry of a tuple at each specific index by splicing the corresponding subtuples
    cd1faces = ([(sx[1:k-1]..., sx[k+1:end]...) for k in eachindex(sx)])
    return cd1faces
end
# This function returns the n-th boundary matrix Dn with coefficients Z/pZ or R of a simplicial complex given as a dimension_dict.
# (p should be a prime. For convenience, p = 0 means R.)
function boundary_matrix(dimension_dict, n, p=0)
    # check if p is 0 or prime
    if p != 0 && !isprime(p)
        throw(DomainError(p, "p should be 0 or a prime."))
    end
    # cases n = 0 and n = dim+1 should be treated separately
    if n == 0
        Dn = zeros(p == 0 ? Float64 : Mod{p}, 1, length(dimension_dict[n]))
    elseif n == maximum(keys(dimension_dict)) + 1
        Dn = zeros(p == 0 ? Float64 : Mod{p}, length(dimension_dict[n-1]), 1)
    else
        Dn = zeros(p == 0 ? Float64 : Mod{p}, length(dimension_dict[n-1]), length(dimension_dict[n]))
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
# determine the mod p Betti numbers of a simplicial complex given via its maximal simplices
function betti_numbers(dimension_dict, p)
    # build the 'dimension dictionary' of a simplicial complex given via its maximal simplices
    #dimension_dict = dim_dict_scx(max_scx)
    # determine the dimension of max_scx
    d = maximum(keys(dimension_dict))
    # initialize the vector of Betti numbers
    betti = []
    # the current (first) boundary matrix
    Dk1 = boundary_matrix(dimension_dict, 0, p)
    # determine the Betti numbers from two consecutive boundary matrices
    for k in 1:(d)
        # the next boundary matrix
        Dk = boundary_matrix(dimension_dict, k, p)
        # perform simultaneous reduction on D_{k-1} and D_k 
        # and store the corresponding Betti number
        push!(betti, simultaneous_reduce(Dk1, Dk))
        # use Dk as the current (first) boundary matrix
        Dk1 = Dk
    end
    return betti
end
function simultaneous_reduce(A, B)
    # throw a 'wrong size' error
    if size(A)[2] != size(B)[1]
        throw(DimensionMismatch("cols(A) != rows(B)"))
    end
    # store the number of rows and columns of A 
    n_rows, n_cols = size(A)
    # perform a column reduction on A and do the inverse operations on rows of B
    i, j = 1, 1
    last_nonzero_col = 0
    while i <= n_rows && j <= n_cols
        # check if A[i, j] can be used as a column pivot ...
        if A[i, j] == 0
            # if not, find a column pivot in the current row
            nonzero_col = j + 1
            while nonzero_col <= n_cols && A[i, nonzero_col] == 0
                nonzero_col += 1
            end
            # if there are only zeros in the current row to the right
            if nonzero_col == n_cols + 1
                # go to a new row
                i += 1
                continue
            end
            # otherwise swap the columns of A ...
            A[:, [j nonzero_col]] = A[:, [nonzero_col j]]
            # ... and corresponding rows of B
            B[[j nonzero_col], :] = B[[nonzero_col j], :]
        end
        # store last nonzero column (for additional row reduction on B later)
        last_nonzero_col = j
        # ... and set A[i, j] as the new column pivot entry 
        pivot = A[i, j]
        # set that column pivot to 1 (in A)
        A[:, j] = A[:, j] / pivot
        # and do the inverse operation on B
        B[j, :] = B[j, :] * pivot
        # annihilate the nonzero elements to the right of the pivot in A
        # and do the opposite transformations on B
        for col in (j+1):n_cols
            scale = A[i, col]
            A[:, col] -= scale * A[:, j]
            B[j, :] += scale * B[col, :]
        end
        # increase row and column indices and search for the next pivot
        i += 1
        j += 1
    end
    # now perform a row reduction on nonzero rows of B 
    # (Performing inverse operations on columns of A would only be required to keep track of the generators...)
    n_rows, n_cols = size(B)
    i, j = last_nonzero_col + 1, 1
    last_nonzero_row = last_nonzero_col
    while i <= n_rows && j <= n_cols
        # check if B[i, j] can be used as a column pivot ...
        if B[i, j] == 0
            # if not, find a row pivot in the current column
            nonzero_row = i + 1
            while nonzero_row <= n_rows && B[nonzero_row, j] == 0
                nonzero_row += 1
            end
            # if there are only zeros in the current column downwards
            if nonzero_row == n_rows + 1
                # go to a new column
                j += 1
                continue
            end
            # otherwise swap the rows of B
            B[[i nonzero_row], :] = B[[nonzero_row i], :]
            # (The corresponding columns of the column-reduced A are 0, no need to swap those...)
        end
        # store last nonzero row (to return the 'Betti number')
        last_nonzero_row = i
        # Set B[i, j] as the new pivot entry ...
        pivot = B[i, j]
        # ... and set that pivot to 1.
        B[i, :] = B[i, :] / pivot
        # annihilate the nonzero elements below the pivot in B
        # (No need for inverse operations on the columns of the reduced A...)
        for row in (i+1):n_rows
            scale = B[row, j]
            B[row, :] -= scale * B[i, :]
        end
        # increase row and column indices and search for the next pivot
        i += 1
        j += 1
    end
    # return the corresponding Betti number
    return n_rows - last_nonzero_row
end
function candidate_radii(MB_edges, MB_triangles, MB_tetrahedra)
    Rs = Float64[]
    append!(Rs, values(MB_edges))
    append!(Rs, values(MB_triangles))
    append!(Rs, values(MB_tetrahedra))
    sort!(unique!(Rs))
end
function covers_sphere(simplex)
    β = betti_numbers(simplex, 2)
    return β[1] == 1 && β[2] == 0 && β[3] == 1
end
function find_R(pts, MB_edges, MB_triangles, MB_tetrahedra)
    Rs = candidate_radii(MB_edges, MB_triangles, MB_tetrahedra)
    for R in Rs
        #println("Trying candidate R = ", R)
        #flush(stdout)  # immediate output
        simplex = cech(R, pts, MB_edges, MB_triangles, MB_tetrahedra)
        #println("Constructed Čech complex for R = ", R)
        #flush(stdout)
        if covers_sphere(simplex)
            return R
        end
    end
    return nothing
end

function plot_edges_triangles(pts, edges, triangles)
    fig = Figure(size=(800, 800))
    ax = Axis3(fig[1, 1], aspect=:data)

    # Sphere mesh
    θ = range(0, 2π, length=80)
    φ = range(0, π, length=40)
    x = [sin(phi) * cos(th) for phi in φ, th in θ]
    y = [sin(phi) * sin(th) for phi in φ, th in θ]
    z = [cos(phi) for phi in φ, th in θ]

    surface!(ax, x, y, z, color=RGBA(0.5, 0.7, 1.0, 0.2), transparency=true)

    # Points
    scatter!(ax, pts[:, 1], pts[:, 2], pts[:, 3], color=RGBA(1, 0, 0, 1), markersize=15)

    # Edges
    for (i, j) in edges
        lines!(ax, [pts[i, 1], pts[j, 1]], [pts[i, 2], pts[j, 2]], [pts[i, 3], pts[j, 3]], color=:black, linewidth=2)
    end

    # Triangles (filled)
    for (i, j, k) in triangles
        poly!(ax, Point3f[pts[i, :], pts[j, :], pts[k, :],], color=RGBA(1.0, 0.6, 0.2, 0.6), transparency=true)
    end

    display(fig)
end


##################
#
#REPORT RENDERING
#
#
#################


function plot_save(pts::Matrix{Float64}, edges, triangles; filename="output.png")
    # Convert points to Float32 tuples
    pts3f = [(Float32(pts[i,1]), Float32(pts[i,2]), Float32(pts[i,3])) for i in 1:size(pts,1)]

    # Create a Scene directly, not via Figure/Axis3
    scene = Scene(resolution = (800, 800), camera = campixel!)

    # Add points
    scatter!(scene, pts3f, color=:red, markersize=10)

    # Add edges
    for (i,j) in edges
        lines!(scene, [pts3f[i], pts3f[j]], color=:black, linewidth=2)
    end

    # Add triangles (filled)
    for (i,j,k) in triangles
        poly!(scene, [pts3f[i], pts3f[j], pts3f[k]], color=RGBA(1.0f0, 0.5f0, 0.0f0, 0.8f0))
    end

    # Set camera position directly
    cam3d!(scene,
           eyeposition = (0f0, 0f0, 0f0),  # negative Z, slightly +Y
           lookat      = (0f0, 0f0, 0f0),
           up          = (0f0, 1f0, 0f0))

    # Force render before saving
    display(scene)  # ensures the camera is applied

    # Save the scene exactly from this angle
    save(filename, scene)
end


##################
#
# This is for Čech coverage with dataset
#
##################

function R_preset()
    pts = read_points("A_sensors_data.txt")
    MB_edges = MiniBall_edges(pts)
    MB_triangles = MiniBall_triangles(pts)
    MB_tetrahedra = MiniBall_tetrahedra(pts)
    R = find_R(pts, MB_edges, MB_triangles, MB_tetrahedra)
    complex = cech(R, pts, MB_edges, MB_triangles, MB_tetrahedra)
    E = complex[1]
    T = complex[2]
    plot_edges_triangles(pts, E, T)
    return R
end


# FINDING OBSOLETE Points
function filter_complex_by_vertices(complex, excluded::Set{Any})
    new_complex = Dict{Int,Vector{Tuple}}()
    for (k, simplices) in complex
        # Keep simplices that do NOT contain any excluded vertex
        new_complex[k] = [c for c in simplices if all(v -> !(v in excluded), c)]
    end
    return new_complex
end
function biggest_excluded(complex, excluded)
    N = length(complex[0])
    E_set = Vector{Vector{Int}}()
    biggest_excl = isempty(excluded) ? 0 : maximum(excluded)
    for n in (biggest_excl+1):N
        println(union(excluded, [n]))
        excluded_n = union(Int.(excluded), [n])
        excluded_set = Set(excluded_n)

        temp_complex = filter_complex_by_vertices(complex, excluded_set)

        β = betti_numbers(temp_complex, 2)
        println(β)
        if β[1] == 1 && β[2] == 0 && β[3] == 1
            deeper = biggest_excluded(temp_complex, excluded_n)
            if deeper !== nothing
                append!(E_set, deeper)
            else
                push!(E_set, union(excluded, [n]))
            end
        end
    end
    return isempty(E_set) ? nothing : E_set
end
function R_extra()
    pts = read_points("A_sensors_data.txt")
    MB_edges = MiniBall_edges(pts)
    MB_triangles = MiniBall_triangles(pts)
    MB_tetrahedra = MiniBall_tetrahedra(pts)
    R = find_R(pts, MB_edges, MB_triangles, MB_tetrahedra)
    #R = 0.4
    complex = cech(R, pts, MB_edges, MB_triangles, MB_tetrahedra)
    be = biggest_excluded(complex, [])


    #plot_edges(pts, E)

    return be
end

function R_excluded(extra_points)
    pts = read_points("A_sensors_data.txt")
    MB_edges = MiniBall_edges(pts)
    MB_triangles = MiniBall_triangles(pts)
    MB_tetrahedra = MiniBall_tetrahedra(pts)
    R = find_R(pts, MB_edges, MB_triangles, MB_tetrahedra)
    println("R done")
    pts2 = pts[setdiff(1:size(pts, 1), extra_points), :]
    MB_edges = MiniBall_edges(pts2)
    MB_triangles = MiniBall_triangles(pts2)
    MB_tetrahedra = MiniBall_tetrahedra(pts2)
    complex = complex = cech(R, pts2, MB_edges, MB_triangles, MB_tetrahedra)
    println(betti_numbers(complex, 2))
    println("rendering")
    E = complex[1]
    T = complex[2]
    plot_edges_triangles(pts2, E, T)
end





# POINT GENERATION
function fibonacci_sphere(N::Integer)
    pts = zeros(Float64, N, 3)
    ϕ = π * (sqrt(5) - 1)             # golden angle in radians
    for k in 0:(N-1)
        z = 1.0 - 2.0 * (k) / ((N))     # from 1 to -1
        r = sqrt(1 - z * z)          # x2 +y2 = 1 - z2, x2+y2=r2
        θ = ϕ * k  # golden angle increment
        x = r * cos(θ)
        y = r * sin(θ)
        pts[k+1, :] = [x, y, z]
    end
    return pts
end

#Thomson optimization
function grad(pts)
    n = size(pts, 1)
    F = zeros(Float64, n, 3)
    for i in 1:(n-1)
        p_i = pts[i, :]
        for j in (i+1):n
            p_j = pts[j, :]
            v = p_i - p_j
            d = norm(p_i - p_j)
            F[i, :] .+= (2 * v / (d^2))
            F[j, :] .-= (2 * v / (d^2))
        end
    end
    return F
end
function relax_on_sphere!(pts; steps=10000, dt=0.005, cooling=0.999, text=false)
    for step in 1:steps
        F = grad(pts)
        pts .+= dt .* F
        for i in 1:size(pts, 1)
            pts[i, :] ./= norm(pts[i, :])
        end
        dt *= cooling

        if text && (step % 50 == 0 || step == 1)
            # measure minimal pairwise distance
            md = minimal_pairwise_distance(pts)
            println("it=$step dt=$(round(dt, sigdigits=3)) min-dist=$(round(md, sigdigits=6))")
        end
    end
end
function minimal_pairwise_distance(pts)
    N = size(pts, 1)
    md = Inf
    for i in 1:N-1
        for j in i+1:N
            d = norm(pts[i, :] - pts[j, :])
            if d < md
                md = d
            end
        end
    end
    return md
end

#HOLE FILLING
function random_point_on_sphere()
    v = randn(3)
    return v / norm(v)
end
function refine_covering!(pts; steps=1000, samples=10000, α=0.02)
    N = size(pts, 1)

    for _ in 1:steps
        worst_x = zeros(3)
        worst_d = -Inf
        worst_i = 1

        for _ in 1:samples
            x = random_point_on_sphere()
            d, i = findmin(norm(x .- pts[j, :]) for j in 1:N)
            if d > worst_d
                worst_d = d
                worst_x .= x
                worst_i = i
            end
        end

        v = pts[worst_i, :] - worst_x
        pts[worst_i, :] -= α * v
        pts[worst_i, :] /= norm(pts[worst_i, :])
    end
end



#######
#
# fibonacci lattice 
#
#######
function r_fib()
    pts = fibonacci_sphere(50)
    r = find_r(pts)
    D = distance_table(pts)
    #V = collect(1:(size(pts,1)))
    E = VR(r, D)
    plot_edges(pts, E)

    return r
end
function R_fib()
    pts = fibonacci_sphere(50)
    MB_edges = MiniBall_edges(pts)
    MB_triangles = MiniBall_triangles(pts)
    MB_tetrahedra = MiniBall_tetrahedra(pts)
    R = find_R(pts, MB_edges, MB_triangles, MB_tetrahedra)
    complex = cech(R, pts, MB_edges, MB_triangles, MB_tetrahedra)
    E = complex[1]
    T = complex[2]
    plot_edges_triangles(pts, E, T)
    return R
end


#######
#
# Thomson
#
#######
function r_thomson()
    pts = fibonacci_sphere(50)
    relax_on_sphere!(pts)
    r = find_r(pts)
    D = distance_table(pts)
    #V = collect(1:(size(pts,1)))
    E = VR(r, D)
    plot_edges(pts, E)

    return r
end
function R_thomson()
    pts = fibonacci_sphere(50)
    relax_on_sphere!(pts)
    MB_edges = MiniBall_edges(pts)
    MB_triangles = MiniBall_triangles(pts)
    MB_tetrahedra = MiniBall_tetrahedra(pts)
    R = find_R(pts, MB_edges, MB_triangles, MB_tetrahedra)
    complex = cech(R, pts, MB_edges, MB_triangles, MB_tetrahedra)
    E = complex[1]
    T = complex[2]
    plot_edges_triangles(pts, E, T)
    return R
end




#######
#
# Hole filling
#
#######
function r_hole()
    r_star = 0
    repetitions = 1
    for i in 1:(repetitions-1)
        println(i)
        pts = fibonacci_sphere(50)
        refine_covering!(pts)
        r = find_r(pts)
        r_star += r
    end
    pts = fibonacci_sphere(50)
    refine_covering!(pts)
    r = find_r(pts)
    r_star += r
    r_star = r_star / repetitions
    D = distance_table(pts)
    #V = collect(1:(size(pts,1)))
    E = VR(r, D)
    plot_edges(pts, E)

    return (r_star, r)
end
function R_hole()
    R_star = 0
    repetitions = 1
    for i in 1:(repetitions-1)
        println(i)
        pts = fibonacci_sphere(50)
        refine_covering!(pts)
        MB_edges = MiniBall_edges(pts)
        MB_triangles = MiniBall_triangles(pts)
        MB_tetrahedra = MiniBall_tetrahedra(pts)
        R = find_R(pts, MB_edges, MB_triangles, MB_tetrahedra)
        R_star += R
    end
    pts = fibonacci_sphere(50)
    refine_covering!(pts)
    MB_edges = MiniBall_edges(pts)
    MB_triangles = MiniBall_triangles(pts)
    MB_tetrahedra = MiniBall_tetrahedra(pts)
    R = find_R(pts, MB_edges, MB_triangles, MB_tetrahedra)
    R_star += R
    R_star = R_star / repetitions
    complex = cech(R, pts, MB_edges, MB_triangles, MB_tetrahedra)
    E = complex[1]
    T = complex[2]
    plot_edges_triangles(pts, E, T)
    return (R_star, R)
end


#######
#
# hole filling on preset
#
#######
function r_hole_preset()
    pts = read_points("A_sensors_data.txt")
    refine_covering!(pts)
    r = find_r(pts)
    D = distance_table(pts)
    #V = collect(1:(size(pts,1)))
    E = VR(r, D)
    plot_edges(pts, E)

    return r
end
function R_hole_preset()
    pts = read_points("A_sensors_data.txt")
    refine_covering!(pts)
    MB_edges = MiniBall_edges(pts)
    MB_triangles = MiniBall_triangles(pts)
    MB_tetrahedra = MiniBall_tetrahedra(pts)
    R = find_R(pts, MB_edges, MB_triangles, MB_tetrahedra)
    complex = cech(R, pts, MB_edges, MB_triangles, MB_tetrahedra)
    E = complex[1]
    T = complex[2]
    plot_edges_triangles(pts, E, T)
    return R
end





#################################
#
#END OF SENSORS
#
#Presentation part 1
#
#################################

function plot_sphere_with_points(ptsM, r)
    pts = [vec(ptsM[i, :]) for i in 1:size(ptsM, 1)]
    fig = Figure(size=(800, 800))
    ax = Axis3(fig[1,1], aspect=:data)

    # Sphere mesh
    θ = range(0, 2π, length=80)
    φ = range(0, π, length=40)
    x = [sin(phi)*cos(th) for phi in φ, th in θ]
    y = [sin(phi)*sin(th) for phi in φ, th in θ]
    z = [cos(phi) for phi in φ, th in θ]

    # Make sphere transparent
    surface!(ax, x, y, z, color=RGBA(0.5, 0.7, 1.0, 1.0), transparency=false)

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

println("starting normal")
pts1=read_points("A_sensors_data.txt")
r1 = find_r(pts1)
E1 = VR(r1, distance_table(pts1))

MB_edges1 = MiniBall_edges(pts1)
MB_triangles1 = MiniBall_triangles(pts1)
MB_tetrahedra1 = MiniBall_tetrahedra(pts1)
R1 = find_R(pts1, MB_edges1, MB_triangles1, MB_tetrahedra1)
complex1 = cech(R1, pts1, MB_edges1, MB_triangles1, MB_tetrahedra1)
E1_ = complex1[1]
T1 = complex1[2]

println("starting fibonacci")
pts2 = fibonacci_sphere(50)
r2 = find_r(pts2)
E2 = VR(r2, distance_table(pts2))

MB_edges2 = MiniBall_edges(pts2)
MB_triangles2 = MiniBall_triangles(pts2)
MB_tetrahedra2 = MiniBall_tetrahedra(pts2)
R2 = find_R(pts2, MB_edges2, MB_triangles2, MB_tetrahedra2)
complex2 = cech(R2, pts2, MB_edges2, MB_triangles2, MB_tetrahedra2)
E2_ = complex2[1]
T2 = complex2[2]

println("starting thomson")
pts3 = fibonacci_sphere(50)
relax_on_sphere!(pts3)
println("points done")
r3 = find_r(pts3)
E3 = VR(r3, distance_table(pts3))

MB_edges3 = MiniBall_edges(pts3)
MB_triangles3 = MiniBall_triangles(pts3)
MB_tetrahedra3 = MiniBall_tetrahedra(pts3)
R3 = find_R(pts3, MB_edges3, MB_triangles3, MB_tetrahedra3)
complex3 = cech(R3, pts3, MB_edges3, MB_triangles3, MB_tetrahedra3)
E3_ = complex3[1]
T3 = complex3[2]

println("starting random")
pts4 = fibonacci_sphere(50)
refine_covering!(pts4)
println("points done")
r4 = find_r(pts4)
E4 = VR(r4, distance_table(pts4))

MB_edges4 = MiniBall_edges(pts4)
MB_triangles4 = MiniBall_triangles(pts4)
MB_tetrahedra4 = MiniBall_tetrahedra(pts4)
R4 = find_R(pts4, MB_edges4, MB_triangles4, MB_tetrahedra4)
complex4 = cech(R4, pts4, MB_edges4, MB_triangles4, MB_tetrahedra4)
E4_ = complex4[1]
T4 = complex4[2]



complex5 = cech(R4, pts1, MB_edges1, MB_triangles1, MB_tetrahedra1)
E5_ = complex5[1]
T5 = complex5[2]



function plot1()#plot points
    plot_edges(pts1, [])
end
function plot2()#plot network
    println("r = ", r1)
    plot_edges(pts1, E1)
end
function plot3()#plot complex
    println("R = ", R1)
    plot_edges_triangles(pts1, E1_, T1)
end

function plot4()#plot covering
    plot_sphere_with_points(pts1, R1)
end

function plot5()#plot covering with R4 a slightly smaller R
    plot_sphere_with_points(pts1, R4)
end

function plot6()#plot complex a bit smaller
    println("complex with a smaller R")
    plot_edges_triangles(pts1, E5_, T5)
end

function plot7() #plot fibbonaci points and edges
    println("Fibbonacci lattice r = ", r2)
    plot_edges(pts2, E2)
end
function plot8() #plot fibbonaci complex
    println("Fibbonacci lattice R = ", R2)
    plot_edges_triangles(pts2, E2_, T2)
end
#plot_edges(pts3, E3)
function plot9() #plot Thomson complex
    println("Thomson simulation r = ", r3)
    println("Thomson simulation R = ", R3)
    plot_edges_triangles(pts3, E3_, T3)
end
#plot_edges(pts4, E4)
function plot10() #plot hole filling complex
    println("Randomized hole filling r = ", r4)
    println("Randomized hole filling R = ", R4)
    plot_edges_triangles(pts4, E4_, T4)
end


