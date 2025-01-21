#%%
include("../julia/Lie_algebra.jl")

#%%
function test_describe()
    definite_algebra = Definite(3)  # SO(3), for example
    indefinite_algebra = Indefinite(1, 3)  # SO(1,3), Lorentz group

    # Describe the algebras
    describe(definite_algebra)  # Output: Definite Lie Algebra: Dimension 3
    describe(indefinite_algebra)  # Output: Indefinite Lie Algebra: Positive 1, Negative 3
    return 0
end
test_describe()
#%% Example: ortho_normal_metric_matrix
#%%
function test_metric(r = 3, s = 2)
    # Example usage
    # r: Number of +1 components
    # s: Number of -1 components
    space = Indefinite(r, s)
    g = ortho_normal_metric_matrix(space)

    # Print the result
    println("The metric tensor matrix g[i, j] is:")
    println(g)

    return 0
end
test_metric()
test_metric(3, 1)
#%% Test function for so(r, s) generator
function test_so_rs_generator(space::Indefinite)
    println(space)
    # Compute the total dimension of the space
    n = dim(space)
    # Iterate over all pairs (mu, nu) with 1 ≤ mu, nu ≤ n
    for mu in 1:n
        for nu in 1:n
            # Call the generator function
            M = so_rs_generator(space, mu, nu)

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
space = Indefinite(2, 1)  # r = 2, s = 1
test_so_rs_generator(space)  # Test with the given space
space = Definite(3)
generate_all_so_generators(space)

#%%
function test_so_generator()
    # Example usage
    vars = Dict(
        "SO3" => (3, 1, 2),
        "SO5" => (5, 1, 2),
        "error" => (5, 1, 0)  # Invalid indices to test error handling
    )

    for label in keys(vars)
        n, i, j = vars[label]
        println("Testing $label with n=$n, i=$i, j=$j")
        try
            space = Definite(n)  # SO(n)
            generator = so_generator(space, i, j)  # Generate the generator
            println("Result:\n$generator")
        catch e
            println("Error for $label: $e")
        end
        println("-"^30)  # Separator for readability
    end

    return 0
end

test_so_generator()
#%%
function test_boost_generator()
    # Define an indefinite space
    space = Indefinite(3, 2)  # r = 3, s = 2, dim = 5
    # Generate a boost generator
    generator = boost_generator(space, 2, 4)
    # Print the result
    println("Boost generator for indices (2, 4):")
    display(generator)
    try
        display(boost_generator(space, 2, 3))  # error
    catch e
        println("Caught error: $e")
    end
    return 0
end
test_boost_generator()
#%%
function show_so_all_commutator(space::Definite)
    n = dim(space)  # Get the dimension of the space
    N = 1:n
    println(space)
    println("Show all commutators of SO(n)")
    for (i, j, k, l) in Iterators.product(N, N, N, N)
        if i == j || k == l
            continue
        else
            # Compute the commutator of two SO(n) generators
            _ans = commutator(so_generator(space, i, j), so_generator(space, k, l))
            
            # Show the indices and the result
            @show (i, j), (k, l)
            display(_ans)
        end
    end

    return 0
end
space = Definite(3)
show_so_all_commutator(space)
#%%
# using Combinatorics
function check_commutative_basis(n::Int = 4)
    println("check_commutative_basis")
    # dimension n such as 4 
    SO_n = Definite(n)
    four_numbers = combinations(1:n, 4)
    # display(four_numbers)
    for num in four_numbers
        for perm in permutations(num)
            print(perm); print(": ");
            product_matrix = (so_generator(SO_n, perm[1], perm[2]) * so_generator(SO_n, perm[3], perm[4]))
            @show isequal(product_matrix, zeros(n, n))
        end
    end
    return 0
end
# four_numbers = combinations(1:4, 4)
# check_commutative_basis(4)
#%%
function test_so4()
    println("Example: Generate and print all generators for SO(4)")
    # n = 4
    SO4 = Definite(4)
    so_4_generators = generate_all_so_generators(SO4)
    display(so_4_generators)
    plus_gen, minus_gen = so_4_lie_algebra(;return_zero = false)

    check_commutative_basis()
    show_products_with(plus_gen, minus_gen)
    return 0
end
test_so4()
so_4_lie_algebra(return_zero = true)
#%%
so_2_1_lie_algebra()
killing_form(structure_so21, 3)
#%%
so_3_1_lie_algebra()

#%%
using Combinatorics
using IterTools

function check_commutator(space::Indefinite)
    println("check_commutator")
    n = dim(space)
    N = 1:n
    itr = IterTools.product(N, N, N, N)
    g(i, j) = ortho_normal_metric(space, i, j)
    M(i, j) = so_rs_generator(space, i, j)
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
            return 1 # Failed
            
        end
    end
    println("Success")
    return 0
end
for (r, s) in [(2, 1), (3, 1)]
    @time check_commutator(Indefinite(r, s))
    
end
#%%
display(sl_generator(SpecialLinear(2), 1, 1))
#%%
sl_2_R_lie_algebra()

#%%

#%%
