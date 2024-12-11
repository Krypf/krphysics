# a non-parallelizable 2-dimensional manifold:
M = Manifold(2, 'M')
U = M.open_subset('U') ; V = M.open_subset('V')
M.declare_union(U,V)   
c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y), intersection_name='W',
                            restrictions1= x>0, restrictions2= u+v>0)
uv_to_xy = xy_to_uv.inverse()
W = U.intersection(V)
eU = c_xy.frame() ; eV = c_uv.frame()
# Differential form of degree 1 on M
a = M.diff_form(1, name='a') ; a
print(a.parent())
# Setting the components of the 1-form in a consistent way:
a[eU,:] = [-y, x]
a.add_comp_by_continuation(eV, W, c_uv)
print(a.display(eU))
print(a.display(eV))
# The exterior derivative of the 1-form is a 2-form:
da = a.exterior_derivative() ; print(da)
print(da.display(eU))
print(da.display(eV))