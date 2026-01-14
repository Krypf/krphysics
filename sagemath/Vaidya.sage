# 1. Initialize the Manifold
# We use 4 dimensions. 'M' is the manifold name.
M = Manifold(4, 'M', structure='Lorentzian')

# 2. Define Coordinates
# u: retarded time, r: radial distance
# th (theta) and ph (phi): angular coordinates
X.<u,r,th,ph> = M.chart(r'u r:(0,oo) th:(0,pi):\theta ph:(0,2*pi):\phi')

# 3. Define the Mass Function m(u)
# We define 'm' as an abstract function of 'u'
m = function('m')(u)

# 4. Define the Metric Tensor g
g = M.metric()

# The metric components:
# g_uu = -(1 - 2m(u)/r)
g[0,0] = -(1 - 2*m/r)

# g_ur = -1 (and g_ru = -1 due to symmetry)
# Note: For INGOING Vaidya (v), this would be +1. For OUTGOING (u), it is -1.
g[0,1] = -1 

# g_thth = r^2
g[2,2] = r^2

# g_phph = r^2 * sin(theta)^2
g[3,3] = r^2 * sin(th)^2

# 5. Compute Curvature Tensors
print("Computing Ricci Tensor...")
Ric = g.ricci()
R_scalar = g.ricci_scalar().display() # 0 

print("Computing Einstein Tensor...")
# G = g.einstein()

# 6. Display the Results
print("\n--- Metric g ---")
g.display()

print("\n--- Einstein Tensor G (Non-zero components) ---")
# The Vaidya metric is a solution to Einstein's equations for a Null Fluid.
# We expect only the G_uu component to be non-zero.
# G.display_comp(only_nonredundant=True)

print("\n--- Calculation of G_uu ---")
# Let's extract the u,u component specifically to see the relation to mass loss
# Val_Guu = G[0,0]
# show(Val_Guu)