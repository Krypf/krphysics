"""
topology_optimized_v2.jl

Enumerates all topologies on X = {0,...,n-1} efficiently using:

  1. Bitmask representation   — subsets of X as UInt32; set ops are bitwise O(1).
  2. DFS with closure pruning — instead of iterating over all 2^(2^n-2) subsets,
                                we build the topology incrementally via DFS.
                                At each step we try adding the next candidate open
                                set, immediately compute the closure, and prune the
                                branch if the closure forces sets that should have
                                been added earlier (Kamae / canonical-extension trick).
  3. Symmetry reduction       — only canonical representatives (lex-min under S_n)
                                are recorded; orbit sizes give the true total count.

Expected counts (OEIS A000798):
  n=1: 1,  n=2: 4,  n=3: 29,  n=4: 355,  n=5: 6942
"""

using Combinatorics   # for permutations()

# ---------------------------------------------------------------------------
# Bitmask utilities
# ---------------------------------------------------------------------------

full_set(n::Int)::UInt32 = UInt32((1 << n) - 1)

function bitmask_to_set(mask::UInt32, n::Int)::Vector{Int}
    return [i for i in 0:n-1 if (mask >> i) & 1 == 1]
end

function print_topology(coll::Vector{UInt32}, n::Int)
    sets = [bitmask_to_set(m, n) for m in sort(coll)]
    print("{ ")
    for (i, s) in enumerate(sets)
        print(isempty(s) ? "∅" : string(s))
        i < length(sets) && print(", ")
    end
    println(" }")
end

# ---------------------------------------------------------------------------
# Closure under finite intersections and unions
# ---------------------------------------------------------------------------

"""
    compute_closure(open_sets_set, n) -> Set{UInt32}

Given a Set of open sets (bitmasks), return the smallest superset closed under
pairwise ∩ and ∪.  Since X is finite, pairwise closure implies full closure.
"""
function compute_closure(open_sets::Set{UInt32}, n::Int)::Set{UInt32}
    current = copy(open_sets)
    changed = true
    while changed
        changed = false
        arr = collect(current)
        for A in arr, B in arr
            for new_set in (A & B, A | B)
                if !(new_set in current)
                    push!(current, new_set)
                    changed = true
                end
            end
        end
    end
    return current
end

# ---------------------------------------------------------------------------
# Symmetry / canonicality
# ---------------------------------------------------------------------------

"""
    apply_perm(mask, perm) -> UInt32

Permute the elements of a bitmask according to `perm` (1-indexed Vector{Int}).
"""
function apply_perm(mask::UInt32, perm::Vector{Int})::UInt32
    result = UInt32(0)
    for (i, j) in enumerate(perm)
        if (mask >> (i-1)) & 1 == 1
            result |= UInt32(1) << (j-1)
        end
    end
    return result
end

"""
    canonical_form(coll_set, perms) -> Vector{UInt32}

Return the lex-minimum image of the sorted collection under all permutations.
`perms` is the precomputed list of all n! permutations (1-indexed).
"""
function canonical_form(coll_set::Set{UInt32}, perms::Vector{Vector{Int}})::Vector{UInt32}
    best = sort(collect(coll_set))
    for perm in perms
        image = sort([apply_perm(m, perm) for m in coll_set])
        if image < best
            best = image
        end
    end
    return best
end

"""
    is_canonical(coll_set, perms) -> Bool

True iff the sorted collection is already the lex-minimum among all its images.
"""
function is_canonical(coll_set::Set{UInt32}, perms::Vector{Vector{Int}})::Bool
    current = sort(collect(coll_set))
    for perm in perms
        image = sort([apply_perm(m, perm) for m in coll_set])
        if image < current
            return false
        end
    end
    return true
end

"""
    orbit_size(coll_set, perms) -> Int

Number of distinct images under all permutations = |orbit| = n! / |stabilizer|.
"""
function orbit_size(coll_set::Set{UInt32}, perms::Vector{Vector{Int}})::Int
    images = Set{Vector{UInt32}}()
    for perm in perms
        push!(images, sort([apply_perm(m, perm) for m in coll_set]))
    end
    return length(images)
end

# ---------------------------------------------------------------------------
# DFS enumeration with pruning
# ---------------------------------------------------------------------------

"""
    dfs_enumerate!(results, current_open, candidates, idx, n, canonical_only, perms)

Depth-first search over subsets of `candidates[idx:end]`.

`current_open` : Set{UInt32} — the open sets chosen so far (always includes ∅ and X,
                 and is already closed under ∩ and ∪).
`candidates`   : Vector{UInt32} — all proper nonempty subsets of X, sorted.
`idx`          : Int — next candidate index to consider (branch or skip).

Pruning rule (canonical extension):
  When we try adding candidates[idx], we first compute the closure of
  current_open ∪ {candidates[idx]}.  Any new set introduced by the closure
  that is NOT among candidates[idx:end] must already be in current_open;
  otherwise the closure would add a set with a smaller index than idx, meaning
  that set should have been added in an earlier recursive call — this branch
  would duplicate a topology already found.  We prune it.
"""
function dfs_enumerate!(
    results       ::Vector{Set{UInt32}},
    current_open  ::Set{UInt32},
    candidates    ::Vector{UInt32},
    idx           ::Int,
    n             ::Int,
    canonical_only::Bool,
    perms         ::Vector{Vector{Int}},
)
    # Record the current topology (it is valid: closed, contains ∅ and X)
    if !canonical_only || is_canonical(current_open, perms)
        push!(results, copy(current_open))
    end

    # Try extending with each remaining candidate in order
    for i in idx:length(candidates)
        c = candidates[i]
        c in current_open && continue   # already present after a previous closure

        # Compute closure of current_open ∪ {c}
        trial = copy(current_open)
        push!(trial, c)
        closed = compute_closure(trial, n)

        # Pruning: every set added by the closure that was not in current_open
        # must have an index >= i in `candidates`.
        # If any such set has index < i, this branch is not a canonical extension
        # and would duplicate a topology found on a different path.
        new_sets = setdiff(closed, current_open)
        valid = true
        for s in new_sets
            # Find the index of s in candidates (if it is a proper nonempty subset)
            pos = findfirst(==(s), candidates)
            if pos !== nothing && pos < i
                valid = false
                break
            end
        end
        valid || continue

        # Recurse with the closed collection
        dfs_enumerate!(results, closed, candidates, i + 1, n, canonical_only, perms)
    end
end

# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

"""
    enumerate_topologies(n; canonical_only=false) -> Vector{Set{UInt32}}

Return all topologies on {0,...,n-1}.
If canonical_only=true, return one representative per homeomorphism class
(i.e. per isomorphism type of the underlying set).
"""
function enumerate_topologies(n::Int; canonical_only::Bool=false)::Vector{Set{UInt32}}
    X          = full_set(n)
    # Proper nonempty subsets, sorted (this order defines the pruning condition)
    candidates = sort([UInt32(k) for k in 1:X-1])
    # All n! permutations (1-indexed vectors)
    perms      = [collect(Int, p) .+ 1 for p in permutations(0:n-1)]

    results = Vector{Set{UInt32}}()
    # Start from the minimal topology {∅, X}
    initial = Set{UInt32}([UInt32(0), X])
    dfs_enumerate!(results, initial, candidates, 1, n, canonical_only, perms)
    return results
end

# ---------------------------------------------------------------------------
# Main demo
# ---------------------------------------------------------------------------

function main()
    println("=" ^ 60)
    println("Topology enumeration  (DFS + closure pruning + symmetry)")
    println("=" ^ 60)

    for n in 1:5
        print("\nn = $n ... ")
        flush(stdout)

        perms = [collect(Int, p) .+ 1 for p in permutations(0:n-1)]

        t_full = @elapsed topologies = enumerate_topologies(n)
        count_full = length(topologies)

        t_canon = @elapsed canonical = enumerate_topologies(n; canonical_only=true)
        count_canon = length(canonical)

        # Verify orbit sizes recover the full count
        orbit_sum = sum(orbit_size(c, n, ) for c in canonical)

        # Use precomputed perms for orbit_size
        orbit_sum2 = sum(orbit_size(c, perms) for c in canonical)

        println("done")
        println("  Total topologies          : $count_full   ($(round(t_full,  digits=4)) s)")
        println("  Canonical representatives : $count_canon   ($(round(t_canon, digits=4)) s)")
        println("  Orbit size check          : orbit_sum = $orbit_sum2  $(orbit_sum2 == count_full ? "✓" : "✗ MISMATCH")")

        if n <= 3
            println("  All topologies on {0,...,$(n-1)}:")
            for (i, top) in enumerate(topologies)
                print("    [$i] ")
                print_topology(collect(top), n)
            end
        end
    end
end

# helper overload so orbit_size works with just (coll_set, n)
function orbit_size(coll_set::Set{UInt32}, n::Int)::Int
    perms = [collect(Int, p) .+ 1 for p in permutations(0:n-1)]
    return orbit_size(coll_set, perms)
end

main()