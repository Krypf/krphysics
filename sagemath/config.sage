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
                com = basis[i].bracket(basis[j])
                print((i, j), com.display())
                print(f"[{basis[i]}, {basis[j]}]")

# print tensors
def show_two_tensor(tensor, dim, array=True, Mathematica=False, result="result.txt"):
    if Mathematica:
        with open(result, "w") as f:
            f.write("# Open Result")
    for i in range(dim):
        for j in range(dim):
            # if i < j:
            if array:
                component = tensor[i][j]
            else:# matrix
                component = tensor[i,j]
            if component != 0:
                show(i, j)
                show(component.display())
                if Mathematica:
                    with open(result, "a") as f:
                        f.write(f"({i}, {j}) ")
                        f.write(str(component.display()))
                        f.write("\n")
                    # print(mathematica_code(x))
    return 0

def show_four_tensor(Riem, dim):
    for i in range(dim):
        for j in range(dim):
            for k in range(dim):
                for l in range(k + 1, dim):
                    component = Riem[i][j][k][l]
                    if component != 0:
                        show(i, j, k, l)
                        show(component.display())
