#%%

# Function to construct a generator of SO(n)
function so_generator(n, i, j)
    # n: Dimension of the space
    # i, j: Indices (1-based indexing, 1 ≤ i < j ≤ n)
    if i < 1 || j > n# || 
        error("Invalid indices: Ensure 1 ≤ i < j ≤ $n")
    end

    # Create an n x n zero matrix
    J = zeros(n, n)
    if i == j
        return J
    end
    # Define the antisymmetric matrix
    J[i, j] = +1 # plus sign
    J[j, i] = -1 # minus
    return J
end
function boost_generator(n, i, j)
    # n: Dimension of the space
    # i, j: Indices (1-based indexing, 1 ≤ i < j ≤ n)
    if i < 1 || j > n
        error("Invalid indices: Ensure 1 ≤ i < j ≤ $n")
    end

    # Create an n x n zero matrix
    J = zeros(n, n)
    
    # Define the antisymmetric matrix
    J[i, j] = +1
    J[j, i] = +1
    return J
end
# Define a function for the commutator
function commutator(A::AbstractMatrix, B::AbstractMatrix)
    return A * B - B * A
end
#%%

function so_all_commutator(n)
    N = 1:n
    for (i, j, k, l) in Iterators.product(N, N, N, N)
        if i == j || k == l
            continue
        else
            _ans = commutator(so_generator(n, i, j), so_generator(n, k, l))
            @show (i, j), (k, l)
            display(_ans)
        end
    end
    return 0
end

# Generate and store all generators of SO(n) in a dictionary
function generate_all_so_generators(n)
    generators = Dict{Tuple{Int, Int}, Matrix{Float64}}()
    for i in 1:(n-1)  # Loop over i, ensuring i < j
        for j in (i+1):n
            generators[(i, j)] = so_generator(n, i, j)
        end
    end
    return generators
end
#%%
function so_4_lie_algebra(_test = true)
    n = 4
    # plus modes
    J1234_plus = - (so_generator(n, 1, 2) + so_generator(n, 3, 4)) / 2
    J1342_plus = - (so_generator(n, 1, 3) - so_generator(n, 2, 4)) /2
    J1423_plus = - (so_generator(n, 1, 4) + so_generator(n, 2, 3)) /2
    plus_gen = [J1234_plus, J1342_plus, J1423_plus]
    @show commutator(plus_gen[1],plus_gen[2]) == plus_gen[3]
    @show commutator(plus_gen[2],plus_gen[3]) == plus_gen[1]
    @show commutator(plus_gen[3],plus_gen[1]) == plus_gen[2]
    # minus modes
    J1234_minus = (so_generator(n, 1, 2) - so_generator(n, 3, 4)) / 2# 1243
    J1342_minus = (so_generator(n, 1, 3) + so_generator(n, 2, 4)) / 2# 1324
    J1423_minus = (so_generator(n, 1, 4) - so_generator(n, 2, 3)) / 2# 1432
    minus_gen = [J1234_minus, J1342_minus, J1423_minus]
    @show commutator(minus_gen[1],minus_gen[2]) == minus_gen[3]
    @show commutator(minus_gen[2],minus_gen[3]) == minus_gen[1]
    @show commutator(minus_gen[3],minus_gen[1]) == minus_gen[2]
    if _test
        return plus_gen, minus_gen
    end
    return 0
end

using Combinatorics
function check_commutative_basis()
    n = 4 # dimension
    s4_iter = permutations(1:n)
    for perm in s4_iter
        matrix = (so_generator(n, perm[1], perm[2]) * so_generator(n, perm[3], perm[4]))
        @show isequal(matrix, zeros(n, n))
    end
    return 0
end

# Compute all products with indices
function show_products_with(plus_gen, minus_gen)
    # m = 4
    for (i, plus) in enumerate(plus_gen)
        for (j, minus) in enumerate(minus_gen)
            @show (((i, j), commutator(plus, minus)))
        end
    end
    return 0
end

function test_so4()
    # Example: Generate and print all generators for SO(4)
    n = 4
    so_n_generators = generate_all_so_generators(n)
    display(so_n_generators)
    plus_gen, minus_gen = so_4_lie_algebra()

    check_commutative_basis()
    show_products_with(plus_gen, minus_gen)
    return 0
end
# test_so4()
#%%
function so_2_1_lie_algebra()
    n = 3
    K13 = boost_generator(n, 1, 3)
    K23 = boost_generator(n, 2, 3)
    J12 = so_generator(n, 1, 2)
    
    @show commutator(K13, K23)
    @show commutator(J12, K13)
    @show commutator(J12, K23)

    return 0
end
# so_2_1_lie_algebra()
#%%
function structure_so21(i,j,k)
    c123 = -1
    # return the same Killing form as a value of c123 = +1
    if (i,j,k) == (1,2,3)
        return c123
    elseif (i,j,k) == (2,1,3)
        return -c123
    elseif (i,j,k) == (3,1,2)
        return -c123
    elseif (i,j,k) == (1,3,2)
        return +c123
    elseif (i,j,k) == (2,3,1)
        return -c123
    elseif (i,j,k) == (3,2,1)
        return +c123
    else return 0
    end
end
function killing_form(structure, dim)
    println("The Killing form of the structure constants")
    _ans = zeros(dim, dim)
    for j in 1:dim
        for k in 1:dim
            # indices
            _ans[j, k] = sum(structure(j, m, n) * structure(k, n, m) for m in 1:dim, n in 1:dim)
        end
    end
    return _ans
end
# killing_form(structure_so21, 3)
#%%
function ortho_normal_metric(r, s, i, j)
    # from 1 to r: +1 if i == j
    # from r + 1 to s: -1 if i == j
    # otherwise 0
    if i == j
        if 1 <= i <= r
            return 1
        elseif r + 1 <= i <= r + s
            return -1
        end
    end
    return 0
end
function ortho_normal_metric_matrix(r, s)
    # Initialize the metric matrix g
    n = r + s # Total dimension
    g = Matrix{Int}(undef, n, n)

    # Fill the matrix using the ortho_normal_metric function
    for i in 1:n
        for j in 1:n
            g[i, j] = ortho_normal_metric(r, s, i, j)
        end
    end

    return g
end

function test_metric()
    # Example usage
    r = 3  # Number of +1 components
    s = 2  # Number of -1 components
    g = ortho_normal_metric_matrix(r, s)

    # Print the result
    println("The metric tensor g[i, j] is:")
    println(g)

    return 0
end
#%%

function kronecker_delta(i, j)
    return i == j ? 1 : 0
end

function so_rs_generator(r, s, mu, nu)
    # n: Dimension of the space
    # i, j: Indices (1-based indexing, 1 ≤ mu < nu ≤ n)
    n = r + s
    if mu < 1 || nu > n# || 
        error("Invalid indices: Ensure 1 ≤ mu < nu ≤ $n")
    end

    # Create an n x n zero matrix
    M = zeros(n, n)
    if mu == nu
        return M
    end
    # Define the generator matrix
    for i in 1:n
        for j in 1:n
            M[i, j] = ortho_normal_metric(r, s, mu, i) * kronecker_delta(nu, j) - ortho_normal_metric(r, s, nu, i) * kronecker_delta(mu, j)
        end
    end
    return M
end

function test_so_rs_generator(r, s)
    # Compute the dimension of the space
    n = r + s

    # Iterate over all pairs (mu, nu) with 1 <= mu, nu <= n
    for mu in 1:n
        for nu in 1:n
            # Call the generator function
            M = so_rs_generator(r, s, mu, nu)

            # Print the resulting matrix
            println("Generator M for mu = $mu, nu = $nu:")
            println(M)

            # Optionally, you can add assertions to validate properties of M
            @assert size(M) == (n, n) "Matrix size is incorrect"
            # @assert M[nu, mu] == -M[mu, nu] "Matrix is antisymmetric" # Test antisymmetry
        end
    end
end

# Example usage:
# test_so_rs_generator(2, 1)  # Test with r = 2, s = 1

#%%
using Combinatorics
using IterTools

function check_commutator(r, s)
    n = r + s
    N = 1:n
    itr = IterTools.product(N, N, N, N)
    g(i, j) = ortho_normal_metric(r, s, i, j)
    M(i, j) = so_rs_generator(r, s, i, j)
    for (i, j, k, l) in itr
        matrix = commutator(M(i, j), M(k, l))
        will_be_com = g(j, k) * M(i, l) + g(i, l) * M(j, k) - g(j, l) * M(i, k) - g(i, k) * M(j, l)
        if isequal(matrix, will_be_com)
            if !isequal(matrix, zeros(n, n))
                @show (i, j, k, l)
                @show (matrix, will_be_com)
            end
            continue
        else
            @show (matrix, will_be_com)
            return 1
            
        end
    end
    return 0
end
# r, s = 2, 1
# r, s = 3, 1
# @time check_commutator(r, s)
#%%
# The Lie algebra of the Lorentz group SO(3, 1)

using Permutations

function levi_civita(indices::Vector{Int})
    # https://chatgpt.com/share/67867b64-3720-800e-9e4b-37d7b07225ae
    if length(Set(indices)) < length(indices)
        return 0  # Indices are not distinct
    else
        return sign(Permutation(indices))
    end
end

function so_3_1_lie_algebra()
    r, s = 3, 1
    n = r + s
    K1 = boost_generator(n, 1, n)
    K2 = boost_generator(n, 2, n)
    K3 = boost_generator(n, 3, n)
    K = [K1, K2, K3]
    
    J32 = so_generator(n, 3, 2)
    J13 = so_generator(n, 1, 3)
    J21 = so_generator(n, 2, 1)
    J = [J32, J13, J21]

    Jplus = (J + 1im * K) / 2
    Jminus = (J - 1im * K) / 2

    m = n - 1
    count_ = 1
    for x in [J, K, Jplus, Jminus]
        for j in 1:m
            @show count_ , j
            display(1im * x[j])# hermitian
        end
        count_ += 1
    end

    for (i, j) in permutations(1:m, 2)
        @show (i, j)
        # @show commutator(J[i], J[j]) == sum(levi_civita([i, j, k]) * J[k] for k in 1:m)
        # @show commutator(K[i], K[j]) == - sum(levi_civita([i, j, k]) * J[k] for k in 1:m)# minus
        # @show commutator(J[i], K[j]) == sum(levi_civita([i, j, k]) * K[k] for k in 1:m)
        # @show commutator(Jplus[i], Jplus[j]) == sum(levi_civita([i, j, k]) * Jplus[k] for k in 1:m)
        # @show commutator(Jminus[i], Jminus[j]) == sum(levi_civita([i, j, k]) * Jminus[k] for k in 1:m)
    end
    
    @show K, J
    return 0
end
so_3_1_lie_algebra()
#%%
function sl_n_R_lie_algebra(n, i, j)
    # n: Dimension of the space
    # i, j: Indices (1-based indexing, 1 ≤ i < j ≤ n)
    N = 1:n
    # Create an n x n zero matrix
    M = zeros(n, n)
    if !(i ∈ N && j ∈ N)
        error("Invalid indices: Ensure i, j = 1,...,n")
    end
    if i == j
        if !(i ∈ 1:(n-1)) error(1) end
        M[i, i] = 1
        M[i+1, i+1] = -1
    else
        M[i, j] = 1
    end
    return M
end

function sl_2_R_lie_algebra()
    # sl(2; R) Lie algebra 
    n = 2; N = 1:n;
    d = Int(n*n - 1); D = 1:d
    diag1 = sl_n_R_lie_algebra(n, 1, 1)
    delta12 = sl_n_R_lie_algebra(n, 1, 2)
    delta21 = sl_n_R_lie_algebra(n, 2, 1)

    standard_basis = [diag1, delta12, delta21]
    println("commutation relations of the standard basis ")
    for (i, j) in combinations(D, 2)
        @show (i, j)
        display(commutator(standard_basis[i], standard_basis[j]))
    end

    normalization = 2
    e1 = diag1 / normalization
    e2 = (delta12 + delta21) / normalization
    e3 = - (delta12 - delta21) / normalization
    basis = [e1, e2, e3]
    for (i, j) in combinations(D, 2)
        @show (i, j)
        LHS = commutator(basis[i], basis[j])
        display(LHS)
        RHS = sum(structure_so21(i, j, k) * basis[k] for k in D) # not for N
        @show RHS
        @show LHS - RHS
    end
    return 0
end
# sl_2_R_lie_algebra()

#%%
function sl2_R_test()
    n = 2
    diag1 = sl_n_R_lie_algebra(n, 1, 1)
    delta12 = sl_n_R_lie_algebra(n, 1, 2)
    delta21 = sl_n_R_lie_algebra(n, 2, 1)
    commutator(diag1, delta21)
    e1 = diag1 / 2
    e2 = (delta12 + delta21) / 2
    e3 = - (delta12 - delta21) / 2
    basis = [e1, e2, e3]
    @show commutator(basis[1], basis[2])-e3
    @show structure_so21(1,2,3)
    return 0
end
#%%
println("Module: julia/Lie_algebra.jl")