#%%
from sympy.abc import u
from sympy import Matrix, Symbol
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
P
#%%
U[0,1] == (psi1 / (phi1 + phi2)).simplify()
#%%
U[1,1] == (psi2 / (phi1 + phi2)).simplify()

# %%
U[1-1,1-1] == - (phi2 * psi1 / (phi1 + phi2)).simplify()
#%%
U[2-1,1-1] == (-1)**2 * (phi1 * psi2 / (phi1 + phi2)).simplify()
# %%

# %%

# %%
