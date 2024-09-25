import sympy as sp
from sympy import eye, Point, Line, Circle, DotProduct
from utils.result import print_md, print_latex

class EuclideanSpace:
    def __init__(self, dimension: int):
        self.dimension = dimension
        # Create Point(0, 0, ..., 0) using generator expression
        # https://chatgpt.com/share/66f30ee0-9064-800e-a4c1-17ebcc4ed730
        self.origin = Point(*(0 for _ in range(dimension))) 
        
    def coefficient_symbols(self, _symbol: str = 'x'):
        x = _symbol + '1:' + str(self.dimension + 1)
        return sp.symbols(x)
    
    def basis_symbols(self, _symbol: str = 'e'):
        x = _symbol + '1:' + str(self.dimension + 1)
        return sp.symbols(x)
    
    def standard_basis(self):
        n = self.dimension
        # Create an n x n identity matrix
        identity_matrix = eye(n)
        # Convert the identity matrix columns to (ordered) standard basis vectors
        return [identity_matrix.col(j) for j in range(n)]
    

def right_line(point_1, point_2):
    return Line(point_1, point_2)

class Place(EuclideanSpace):
    def __init__(self, dimension, position=None):
        super().__init__(dimension) # Inherit dimensions and origin from EuclideanSpace
        self.position = position
        
    def coefficients_motion(self, function_name: str = "x", args = sp.Symbol("t")):
        from sympy import Function, Matrix
        d = self.dimension
        # Check if args is a string
        if isinstance(args, str):
            args = sp.Symbol(args)
        return Matrix([Function(f"{function_name}_{i}")(args) for i in range(1, d + 1)])
    
    def moving_frame(self, function_name: str = "e", args = sp.Symbol("t")):
        from sympy import Function, Matrix
        d = self.dimension
        # Check if args is a string
        if isinstance(args, str):
            args = sp.Symbol(args)
        return Matrix([[Function(f"{function_name}_{i}")(args) for i in range(1, d + 1)]])
    
    def fictitious_force(self, args = sp.Symbol("t")):
        x = self.coefficients_motion()
        e = self.moving_frame()
        position_vector = (e * x)[0, 0]
        return position_vector.diff(args, 2)
        

class PhysicalQuantity:
    def __init__(self,
        density=None,
        volume=None,
        velocity=None,
        mass=None,
        acceleration=None,
    ):
        self.density = density
        self.volume = volume
        self.velocity = velocity
        self.mass = mass
        self.acceleration = acceleration
        
    def quantity_matter(self):
        if self.mass:
            return self.mass
        else:        
            return self.density * self.volume
        
    def quantity_motion(self):
        return self.velocity * self.quantity_matter()
    def innate_force(self):
        return 
    
    def impressed_force(self):
        return self.quantity_matter() * self.acceleration
    
    
def law_motion_1st():
    print_md("Every body perseveres in its state of rest, or of uniform motion in a right line, unless it is compelled to change that state by forces impressed thereon.")

def law_motion_2nd():
    from sympy.abc import m, v, F, t
    x = PhysicalQuantity(mass=m, velocity=v)
    p = x.quantity_motion()
    alteration_motion = sp.Function("Delta")(p)
    coefficient_motion = sp.Function("Delta")(t)
    eq_motion = sp.Eq(alteration_motion, F * coefficient_motion)
    # Save the result
    print_md("# LAW II.\nThe alteration of motion is ever proportional to the motive force impressed; and is made in the direction of the right line in which that force is impressed.")
    # Get the relative path of the script
    import os
    script_name = os.path.relpath(__file__) if '__file__' in globals() else 'unknown file'
    print_latex(sp.latex(eq_motion), script_name=script_name)


def law_motion_3rd():
    print("To every action, there is always opposed an equal reaction: or the mutual actions of two bodies upon each other are always equal, and directed to contrary parts.")
    