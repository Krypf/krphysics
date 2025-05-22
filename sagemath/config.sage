def show_vectors(vectors):
    print("Show vectors")
    for v in vectors:
        print(v.display())
        
def show_diffs(forms):
    print("Exterior derivatives of forms")
    for v in forms:
        print(diff(v).display())

def pairing_forms_vectors(forms, vectors):
    for i in range(len(forms)):  # Iterate over indices of forms
        for j in range(len(vectors)):  # Iterate over indices of vectors
            print((i, j), forms[i](vectors[j]).display())  # Apply forms[i] to vectors[j] and display the result

def show_commutators(basis):
    n = len(basis)
    print("Commutators of a basis.")
    for i in range(n):
        for j in range(n):
            if i < j:
                x = basis[i].bracket(basis[j])
                print((i, j), x.display())
                print(f"[{basis[i]}, {basis[j]}]")
