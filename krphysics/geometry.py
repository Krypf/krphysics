"""
Topology T of a set X has
1. the empty set {} and the whole set X
2. The intersection of A in T and B and T belongs to T
3. Any union of elements of T belongs to T
"""
# https://chatgpt.com/share/67027011-d11c-800e-8af9-2f62adcbe941
class Topology:
    """
    Initialize a topology with:
    space: The set of elements (as a Python set)
    collection: The collection of subsets (as a list of sets) satisfying topology axioms
    """
    def __init__(self, space, collection):
        self.space = space
        self.collection = collection
    def has_empty(self):
        return set() in self.collection
    def has_whole(self):
        return self.space in self.collection

    def closure_intersection(self):
    # Axiom 2: collection are closed under intersection (for finite intersections)
        for A in self.collection:
            for B in self.collection:
                if A.intersection(B) not in self.collection:
                    return False
        return True
    
    def closure_all_union(self):
        # Axiom 3: collection are closed under arbitrary union
        for i in range(1, len(self.collection) + 1):
            for union_subset in self._get_all_unions(i):
                if union_subset not in self.collection:
                    return False
        return True

    def _get_all_unions(self, n):
        """ Helper method to compute unions of collection of subsets of length n """
        from itertools import combinations
        return [set.union(*combo) for combo in combinations(self.collection, n)]

    def is_topology(self):
        """ Check if collection satisfy the axioms of a topology """
        if not (self.has_empty() and self.has_whole()):
            print("Warning: The topology must include the empty set and the whole set.")
            return False
    
        if not self.closure_intersection():
            print("Warning: The topology must be closed under intersection.")
            return False
        
        if not self.closure_all_union():
            print("Warning: The topology must be closed under arbitrary union.")
            return False
        
        return True

if __name__ == '__main__':
    # Example usage
    X = {1, 2, 3}
    T = [set(), {1}, {2}, {1, 2}, X]

    top = Topology(X, T)
    print("Is this a valid topology?", top.is_topology())
    
