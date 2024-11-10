#%%
# from sympy.abc import u
from sympy import Matrix, Symbol
from utils.result import Result
result = Result(
    script_name="BF00757025/section3.py",
    file_name="section3.md",
    dir_path="BF00757025"
)
result.create_result()

#%% 
Uinv = Matrix([
    [Symbol("{}^1u_1"), Symbol("{}^1u_2")],
    [Symbol("{}^2u_1"), Symbol("{}^2u_2")],
])
U = Uinv**-1
#%%
psi1 = -1 / Symbol("{}^1u_1")
psi2 = 1 / Symbol("{}^1u_2")
phi1 = - Symbol("{}^2u_1") / Symbol("{}^1u_1")
phi2 = Symbol("{}^2u_2") / Symbol("{}^1u_2")

P = Matrix([
    [psi1,  psi2],
    [phi1,  phi2],
])
result.print_md("## Equations (18) and (19)")
result.print_md("Let us suppose $P = $", end=" ")
result.print_latex(P)
#%%
def print_both_sides(A, B):
    result.print_latex(A)
    result.print_latex(B)
    result.print_md(f"A == B is {A == B}.")
    
result.print_md("## Equations (21) and (22)")
A = U[0,1]
B = (psi1 / (phi1 + phi2)).simplify()
result.print_md("$g^{11} = {}_2 u^{1} = $", end=" ")
print_both_sides(A, B)
#%%
A = U[1,1]
B = (psi2 / (phi1 + phi2)).simplify()
result.print_md("$g^{22} = {}_2 u^{2} = $", end=" ")
print_both_sides(A, B)
# %%
result.print_md("## Equations (23) and (24)")
a = 1; b = 1;
A = U[b-1,a-1]
B = - (phi2 * psi1 / (phi1 + phi2)).simplify()
result.print_md("$K^{11} = {}_a u^{b} = $", end=" ")
print_both_sides(A, B)
#%%
a = 1; b = 2;
A = U[b-1,a-1] 
B = (-1)**2 * (phi1 * psi2 / (phi1 + phi2)).simplify()
result.print_md("$K^{22} = {}_a u^{b} = $", end=" ")
print_both_sides(A, B)
# %%
# K_{2}^{11}
# %% (24)

# %%
