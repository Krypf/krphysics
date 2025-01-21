# https://doc.sagemath.org/html/en/reference/manifolds/sage/manifolds/differentiable/vectorfield.html
class Sphere:
    dim = 2
    def __init__(self,
        args_name = 'theta phi',
        manifold = Manifold(dim, 'Sphere', structure='Riemannian'),
        radius_name = 'a',
        neighborhood_name = 'U',
        ):
        self.radius = var('a')
        # Define the manifold
        self.manifold = manifold
        self.args_name = args_name
        self.args = var(args_name)# = theta, phi 
        # the open set covered by spherical coord.
        self.neighborhood = manifold.open_subset(neighborhood_name) # = U

XS.<theta, phi> = U.chart(r'θ:(0,pi):\theta ϕ:(0,2*pi):\phi')
eth, eph = XS.frame()
wth = S2.one_form(1, 0, name='wth')
wph = S2.one_form(0, sin(theta), name='wth')
e1 = eth
e2 = eph / sin(theta)
commutator_e1_e2 = e1.bracket(e2)
print(commutator_e1_e2.display())