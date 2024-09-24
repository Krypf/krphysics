#%%
from krphysics.classical_mechanics import law_motion_2nd, EuclideanSpace
#%%

if __name__ == '__main__':
    print("classical mechanics")
    d = 3
    E = EuclideanSpace(d)
    x = E.coefficient_symbols('x')
    print(x)
    # law_motion_1st()
    # law_motion_2nd()

#%%
