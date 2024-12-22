# Function to construct a generator of SO(n)
function so_generator(n, i, j)
    # n: Dimension of the space
    # i, j: Indices (1-based indexing, 1 ≤ i < j ≤ n)
    if i >= j || i < 1 || j > n
        error("Invalid indices: Ensure 1 ≤ i < j ≤ $n")
    end

    # Create an n x n zero matrix
    J = zeros(n, n)
    
    # Define the antisymmetric matrix
    J[i, j] = -1 # minus sign
    J[j, i] = +1 # plus
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
# Define a function for the commutator
function commutator(A::AbstractMatrix, B::AbstractMatrix)
    return A * B - B * A
end
#%%

# Example: Generate and print all generators for SO(4)
n = 4
generators = generate_all_so_generators(n)
for ((i, j), gen) in generators
    println("Generator ($i, $j):")
    display(gen)
end
#%%
J1234_plus = (so_generator(n, 1, 2) + so_generator(n, 3, 4)) / 2
J1342_plus = (so_generator(n, 1, 3) - so_generator(n, 2, 4)) /2
J1423_plus = (so_generator(n, 1, 4) + so_generator(n, 2, 3)) /2
plus_gen = [J1234_plus, J1342_plus, J1423_plus]
commutator(plus_gen[1],plus_gen[2]) == plus_gen[3]
commutator(plus_gen[2],plus_gen[3]) == plus_gen[1]
commutator(plus_gen[3],plus_gen[1]) == plus_gen[2]
#%%
J1234_minus = (so_generator(n, 1, 2) - so_generator(n, 3, 4)) / 2# 1243
J1342_minus = (so_generator(n, 1, 3) + so_generator(n, 2, 4)) / 2# 1324
J1423_minus = (so_generator(n, 1, 4) - so_generator(n, 2, 3)) / 2# 1432
minus_gen = [-J1234_minus, -J1342_minus, -J1423_minus]
commutator(minus_gen[1],minus_gen[2]) == minus_gen[3]
commutator(minus_gen[2],minus_gen[3]) == minus_gen[1]
commutator(minus_gen[3],minus_gen[1]) == minus_gen[2]

#%% commutative basis
(so_generator(n, 1, 2) * so_generator(n, 3, 4))
commutator(so_generator(n, 1, 2), so_generator(n, 3, 4))
#%%
# Compute all products with indices
function generate_products_with_indices(plus_gen, minus_gen)
    products = []
    for (i, plus) in enumerate(plus_gen)
        for (j, minus) in enumerate(minus_gen)
            @show (((i, j), commutator(plus, minus)))
        end
    end
    return nothing
end
generate_products_with_indices(plus_gen, minus_gen)
#%%
n = 3
K13 = boost_generator(n, 1, 3)
K23 = boost_generator(n, 2, 3)
J12 = so_generator(n, 1, 2)

@show commutator(K13, K23)
@show commutator(J12, K13)
@show commutator(J12, K23)
#%%
J23 = so_generator(n, 2, 3)
J13 = so_generator(n, 1, 3)
@show commutator(J12, J23)
