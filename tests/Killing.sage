# 1. Define a 3-dimensional differentiable manifold
M = Manifold(3, 'M')

# 2. Define the Cartesian coordinate chart (x, y, z)
C_cart.<x,y,z> = M.chart()

# 3. Define the Spherical coordinate chart (r, th, ph)
# th = theta (polar angle), ph = phi (azimuthal angle)
C_sph.<r,th,ph> = M.chart(r'r:(0,+oo) th:(0,pi):periodic ph:(0,2*pi):periodic')

# 4. Define the transition map from Spherical to Cartesian
# x = r*sin(th)*cos(ph), y = r*sin(th)*sin(ph), z = r*cos(th)
C_sph.transition_map(C_cart, [r*sin(th)*cos(ph), 
                              r*sin(th)*sin(ph), 
                              r*cos(th)])

# Get the basis vectors of the Cartesian chart
e_cart = C_cart.frame()

# 5. Define the Translation Killing Vectors (Dx, Dy, Dz)
# These represent partial derivatives with respect to x, y, and z
Dx = M.vector_field(1, 0, 0, basis=e_cart, name='Dx')
Dy = M.vector_field(0, 1, 0, basis=e_cart, name='Dy')
Dz = M.vector_field(0, 0, 1, basis=e_cart, name='Dz')

# 6. Define the Rotation Killing Vectors (Lx, Ly, Lz)
# Lz = x*dy - y*dx
Lz = x*Dy - y*Dx
Lz.set_name('Lz')

# Lx = y*dz - z*dy
Lx = y*Dz - z*Dy
Lx.set_name('Lx')

# Ly = z*dx - x*dz
Ly = z*Dx - x*Dz
Ly.set_name('Ly')

# 7. Display the results in Spherical Coordinates
print("--- Translation Vectors in Spherical Coordinates ---")
print((Dx.display(C_sph)))
print((Dy.display(C_sph)))
print((Dz.display(C_sph)))

print("\n--- Rotation Vectors in Spherical Coordinates ---")
print((Lx.display(C_sph)))
print((Ly.display(C_sph)))
print((Lz.display(C_sph)))