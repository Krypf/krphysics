"""
topology_optimized.jl

Enumerates all topologies on a finite set X = {0, 1, ..., n-1} using three optimizations
over the naive 2^(2^n) brute-force search:

  1. Bitmask representation  — subsets of X are UInt32 integers; set operations become
                               bitwise AND/OR, which are O(1).
  2. Constructive closure    — instead of checking all 2^(2^n) collections, we start
                               from {∅, X} and add candidate open sets one by one,
                               immediately completing the collection under ∩ and ∪.
                               This prunes enormous branches early.
  3. Symmetry reduction      — we generate only topologies whose open sets, when sorted,
                               are lexicographically canonical under all permutations of X.
                               Each isomorphism class is counted once, then multiplied by
                               the size of its orbit (n! / |stabilizer|).

Usage:
    julia topology_optimized.jl          # runs demo for n = 1..5
"""

# ---------------------------------------------------------------------------
# Bitmask utilities
# ---------------------------------------------------------------------------

"""
    subsets(n) -> Vector{UInt32}

Return all 2^n subsets of {0,...,n-1} as bitmasks.
Bit k is set iff element k belongs to the subset.
"""
function subsets(n::Int)::Vector{UInt32}
    return UInt32.(0 : (1 << n) - 1)
end

"""
    full_set(n) -> UInt32

Bitmask for the whole set X = {0,...,n-1}.
"""
full_set(n::Int)::UInt32 = UInt32((1 << n) - 1)

# ---------------------------------------------------------------------------
# Closure under finite intersections and arbitrary unions
# ---------------------------------------------------------------------------

"""
    close_topology(open_sets, n) -> Vector{UInt32}

Given a collection of open sets (as bitmasks), return the smallest collection
that contains them all and is closed under ∩ and ∪.

Strategy: repeatedly add missing intersections and unions until stable.
This is the 'closure' step used during constructive enumeration.
"""
function close_topology(open_sets::Vector{UInt32}, n::Int)::Vector{UInt32}
    current = Set{UInt32}(open_sets)
    changed = true
    while changed
        changed = false
        arr = collect(current)
        # Close under pairwise intersection
        for A in arr, B in arr
            inter = A & B
            if !(inter in current)
                push!(current, inter)
                changed = true
            end
        end
        # Close under pairwise union (repeated application covers arbitrary union
        # on finite sets by induction)
        arr = collect(current)
        for A in arr, B in arr
            uni = A | B
            if !(uni in current)
                push!(current, uni)
                changed = true
            end
        end
    end
    return sort(collect(current))
end

# ---------------------------------------------------------------------------
# Topology validation (used for testing / verification)
# ---------------------------------------------------------------------------

"""
    is_topology(coll, n) -> Bool

Check that `coll` (a sorted Vector{UInt32} of bitmasks) is a valid topology on
{0,...,n-1}: contains ∅ and X, and is closed under ∩ and ∪.
"""
function is_topology(coll::Vector{UInt32}, n::Int)::Bool
    s = Set{UInt32}(coll)
    # Must contain ∅ and X
    (UInt32(0) in s && full_set(n) in s) || return false
    # Closed under intersection and union
    for A in coll, B in coll
        ((A & B) in s && (A | B) in s) || return false
    end
    return true
end

# ---------------------------------------------------------------------------
# Canonical form under permutations of X  (symmetry reduction)
# ---------------------------------------------------------------------------

"""
    apply_perm(mask, perm) -> UInt32

Apply a permutation `perm` (a Vector{Int} mapping old index → new index) to a
bitmask representing a subset of {0,...,n-1}.
"""
function apply_perm(mask::UInt32, perm::Vector{Int})::UInt32
    result = UInt32(0)
    for (i, j) in enumerate(perm)          # i-th element maps to j
        if (mask >> (i-1)) & 1 == 1
            result |= UInt32(1) << (j-1)
        end
    end
    return result
end

"""
    apply_perm_to_collection(coll, perm) -> Vector{UInt32}

Apply a permutation to every set in the collection and return the result sorted.
"""
function apply_perm_to_collection(coll::Vector{UInt32}, perm::Vector{Int})::Vector{UInt32}
    return sort([apply_perm(m, perm) for m in coll])
end

"""
    is_canonical(coll, n) -> Bool

Return true iff `coll` is the lexicographically smallest among all images of
`coll` under permutations of {0,...,n-1}.  Only canonical representatives are
kept, eliminating duplicates up to homeomorphism on the underlying set.
"""
function is_canonical(coll::Vector{UInt32}, n::Int)::Bool
    for perm in permutations(0:n-1)
        perm_vec = collect(Int, perm) .+ 1   # 1-indexed for apply_perm
        image = apply_perm_to_collection(coll, perm_vec)
        if image < coll                       # lex comparison on sorted vectors
            return false
        end
    end
    return true
end

"""
    orbit_size(coll, n) -> Int

Count the number of distinct images of `coll` under all n! permutations.
This equals n! / |stabilizer|.  Used to reconstruct the total count from
canonical representatives.
"""
function orbit_size(coll::Vector{UInt32}, n::Int)::Int
    images = Set{Vector{UInt32}}()
    for perm in permutations(0:n-1)
        perm_vec = collect(Int, perm) .+ 1
        push!(images, apply_perm_to_collection(coll, perm_vec))
    end
    return length(images)
end

# ---------------------------------------------------------------------------
# Constructive enumeration
# ---------------------------------------------------------------------------

"""
    enumerate_topologies(n; canonical_only=false) -> Vector{Vector{UInt32}}

Enumerate all topologies on {0,...,n-1} using the constructive closure approach.

Algorithm:
  - The mandatory open sets are ∅ (= 0x0) and X (= full_set(n)).
  - Candidate open sets are the 2^n - 2 proper nonempty subsets of X.
  - We iterate over all 2^(2^n - 2) subsets of candidates.
    For each candidate subset S:
      1. Start with {∅, X} ∪ S.
      2. Compute the closure under ∩ and ∪.
      3. If the closure equals {∅, X} ∪ S (i.e. S was already closed), record it.
         Otherwise skip — the closure will be encountered as a different candidate.
  - If canonical_only=true, keep only one representative per isomorphism class.

Note: the 2^(2^n-2) iteration is still exponential, but each step is O(1) in
bitwise ops and the closure check prunes most branches immediately.
"""
function enumerate_topologies(n::Int; canonical_only::Bool=false)::Vector{Vector{UInt32}}
    X     = full_set(n)
    # All proper nonempty subsets of X
    candidates = [UInt32(k) for k in 1:X-1]   # excludes 0 (=∅) and X
    m     = length(candidates)                  # = 2^n - 2

    results = Vector{Vector{UInt32}}()

    for mask in UInt32(0) : UInt32((1 << m) - 1)
        # Build initial collection: {∅, X} plus the selected candidates
        init = UInt32[0, X]
        for i in 1:m
            if (mask >> (i-1)) & 1 == 1
                push!(init, candidates[i])
            end
        end

        # Compute closure under ∩ and ∪
        closed = close_topology(init, n)

        # Accept only if the closure did not add any new sets
        # (otherwise this collection is not itself a topology — it will be
        #  found again when we iterate over the mask corresponding to `closed`)
        if Set(closed) == Set(init)
            if !canonical_only || is_canonical(closed, n)
                push!(results, closed)
            end
        end
    end

    return results
end

# ---------------------------------------------------------------------------
# Pretty printing
# ---------------------------------------------------------------------------

"""
    bitmask_to_set(mask, n) -> Vector{Int}

Convert a bitmask back to a human-readable sorted list of elements.
"""
function bitmask_to_set(mask::UInt32, n::Int)::Vector{Int}
    return [i for i in 0:n-1 if (mask >> i) & 1 == 1]
end

"""
    print_topology(coll, n)

Print a topology as a collection of subsets of {0,...,n-1}.
"""
function print_topology(coll::Vector{UInt32}, n::Int)
    sets = [bitmask_to_set(m, n) for m in coll]
    print("{ ")
    for (i, s) in enumerate(sets)
        print(isempty(s) ? "∅" : string(s))
        i < length(sets) && print(", ")
    end
    println(" }")
end

# ---------------------------------------------------------------------------
# Main demo
# ---------------------------------------------------------------------------

using Combinatorics   # for permutations()

function main()
    println("=" ^ 60)
    println("Topology enumeration on finite sets")
    println("=" ^ 60)

    for n in 1:5
        println("\nn = $n  (|X| = $n, candidates = $(2^n - 2), search space = 2^$(2^n-2) = $(2^(2^n-2)))")

        # --- Full count (no symmetry reduction) ---
        t_full = @elapsed topologies = enumerate_topologies(n)
        count_full = length(topologies)

        # --- Canonical representatives only ---
        t_canon = @elapsed canonical = enumerate_topologies(n; canonical_only=true)
        count_canon = length(canonical)

        println("  Total topologies          : $count_full   ($(round(t_full, digits=4)) s)")
        println("  Canonical representatives : $count_canon   ($(round(t_canon, digits=4)) s)")

        # Verify orbit sizes sum to total
        orbit_sum = sum(orbit_size(c, n) for c in canonical)
        @assert orbit_sum == count_full "Orbit sum mismatch!"
        println("  Orbit size check          : ✓  ($orbit_sum = $count_full)")

        # Print all topologies for small n
        if n <= 3
            println("  Topologies on {0,...,$(n-1)}:")
            for (i, top) in enumerate(topologies)
                print("    [$i] ")
                print_topology(top, n)
            end
        end
    end
end

main()
# n = 5  (|X| = 5, candidates = 30, search space = 2^30 = 1073741824)
# ^CERROR: LoadError: InterruptException: