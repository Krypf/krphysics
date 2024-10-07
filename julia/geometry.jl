module TopologyModule

using Combinatorics  # Import the Combinatorics package

export Topology, has_empty, has_whole, closure_intersection, closure_all_union, get_all_unions, is_topology

struct Topology
    space::Set{Any}
    collection::Set{Set{Any}}
end

function has_empty(t::Topology)
    return Set() in t.collection
end

function has_whole(t::Topology)
    return t.space in t.collection
end

function closure_intersection(t::Topology)
    # Axiom 2: collection is closed under intersection (for finite intersections)
    for A in t.collection
        for B in t.collection
            if !(intersect(A, B) in t.collection)
                return false
            end
        end
    end
    return true
end

function closure_all_union(t::Topology)
    # Axiom 3: collection is closed under arbitrary union
    for i in 1:length(t.collection)
        for union_subset in get_all_unions(t.collection, i)
            if !(union_subset in t.collection)
                return false
            end
        end
    end
    return true
end

function get_all_unions(collection::Set{Set{Any}}, n::Int)
    """ Helper method to compute unions of subsets of length n """
    return [reduce(union, combo) for combo in combinations(collect(collection), n)]
end

function is_topology(t::Topology)
    """ Check if collection satisfies the axioms of a topology """
    if !(has_empty(t) && has_whole(t))
        println("Warning: The topology must include the empty set and the whole set.")
        return false
    end

    if !closure_intersection(t)
        println("Warning: The topology must be closed under intersection.")
        return false
    end

    if !closure_all_union(t)
        println("Warning: The topology must be closed under arbitrary union.")
        return false
    end

    return true
end

end  # module

# Function to convert Vector{Vector{Int64}} to Set{Set{Int64}}
function parse_to_set_of_sets(vectors::Vector{Vector{Int64}})
    return Set(Set(vec) for vec in vectors)
end

# Example usage
function test()
    # X = Set([1, 2])
    # X = Set([1, 2, 3])
    X = Set([1, 2, 3, 4])
    T = Set([Set([]), Set([1]), Set([2]), Set([1, 2]), X])

    top = Topology(X, T)
    println("Is this a valid topology? ", is_topology(top))

    p = powerset(collect(X))
    for t in powerset(collect(p))
        t = parse_to_set_of_sets(t)
        top = Topology(X,t)
        if is_topology((top))
            println(t)
        else
            continue
        end
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    #Test code
    using .TopologyModule
    using Combinatorics
    @time test()
end

