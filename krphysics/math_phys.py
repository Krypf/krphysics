import sympy as sp

class VectorField:
    def __init__(self, components, variables):
        """
        Initialize a vector field.

        Parameters:
        components (list): A list of vector field components [F_x, F_y, F_z].
        variables (list): A list of variables [x, y, z].
        """
        if len(components) != len(variables):
            raise ValueError("The number of components must match the number of variables.")
        
        self.components = components
        self.variables = variables

    def divergence(self):
        """
        Calculate the divergence of the vector field.
        
        Returns:
        sympy.Expr: The divergence of the vector field.
        """
        div = sum(sp.diff(self.components[i], self.variables[i]) for i in range(len(self.variables)))
        return div

    def curl(self):
        """
        Calculate the curl (rotation) of the vector field in 3D.
        
        Returns:
        list: A list representing the components of the curl vector.
        """
        if len(self.variables) != 3:
            raise ValueError("Curl is only defined for 3D vector fields.")
        
        curl_x = sp.diff(self.components[2], self.variables[1]) - sp.diff(self.components[1], self.variables[2])
        curl_y = sp.diff(self.components[0], self.variables[2]) - sp.diff(self.components[2], self.variables[0])
        curl_z = sp.diff(self.components[1], self.variables[0]) - sp.diff(self.components[0], self.variables[1])
        
        return [curl_x, curl_y, curl_z]

    def gradient(self, scalar_field):
        """
        Calculate the gradient of a scalar field.
        
        Parameters:
        scalar_field (sympy.Expr): A scalar function.
        
        Returns:
        list: A list representing the components of the gradient vector.
        """
        grad = [sp.diff(scalar_field, var) for var in self.variables]
        return grad
