###########################################################
#  SphericalMetric.sage
#  A general spherically symmetric spacetime metric class
#  (in the same design as Schwarzschild class)
###########################################################

load("~/krphysics/sagemath/spacetime.sage")

class SphericalMetric(Spacetime):
    def __init__(self, manifold,
                 f_r=None, g_r=None,
                 args_name='t r theta phi',
                 neighborhood_name='U'):
        # 初期化
        super().__init__(manifold)
        self.args_name = args_name
        self.args = var(args_name)  # (t, r, theta, phi)
        self.neighborhood = manifold.open_subset(neighborhood_name)
        
        # デフォルト関数設定
        r = self.args[1]
        if f_r is None:
            self.f_r = function('f')(r)
        elif f_r == "exp":
            self.f_r = exp(2 * function('mu', imag_part_func=0)(r)) # https://ask.sagemath.org/question/61315/assume-a-function-is-real-valued/
        else:
            self.f_r = f_r
            # exit(1)
        if g_r is None:
            self.g_r = function('g')(r, domain='real')
        elif g_r == "exp":
            # self.g_r = exp(2 * function('nu', imag_part_func=0)(r))
            self.g_r = exp(- 2 * function('mu', imag_part_func=0)(r))
            # assume(self.g_r, 'real')
        else:
            self.g_r = g_r
            # exit(1)
        # 符号行列（Lorentz 符号）
        self.sign_matrix = (+1) * diagonal_matrix([-1, 1, 1, 1])

    # 球座標チャート
    def spherical_chart(self):
        return self.neighborhood.chart(self.args_name)

    # 座標フレームと双対フレーム
    def coordinate_frame(self):
        return self.spherical_chart().frame()

    def dual_frame(self):
        return self.spherical_chart().coframe()

    # 成分とポテンシャル関数
    def components(self, number_radius=1, number_theta=2):
        r = self.args[number_radius]
        theta = self.args[number_theta]
        return r, theta, self.f_r, self.g_r

    def potential_functions(self):
        return self.f_r, self.g_r
    
    # 計量テンソル
    def metric_tensor(self):
        r, theta, f_r, g_r = self.components()
        U = self.spherical_chart().domain() # self.neighborhood
        gU = U.metric('g_U')

        gU[0,0] = -f_r
        gU[1,1] = g_r
        gU[2,2] = r^2
        gU[3,3] = r^2 * sin(theta)^2

        return gU
    
    def electromagnetic(self):
        # dom = self.chart.domain()
        x = self.spherical_chart()
        F = x.domain().diff_form(2)
        E = [function('Er')(*x), 0, 0]
        for i in range(1, 3+1):
            F[i,0] = E[i-1] / c
        # B = [function('B1')(*x), function('B2')(*x), function('B3')(*x)]
        F[1,2] = 0 # function('B12')(*x) # B[3-1]  # dx∧dy
        F[2,3] = 0 # function('B23')(*x) # B[1-1]  # dy∧dz
        F[3,1] = 0 # function('B31')(*x) # B[2-1]  # dz∧dx
        return F

    def excitation(self):
        x = self.spherical_chart()
        # dom = self.chart.domain()
        exc = x.domain().diff_form(2)
        g_inv = self.metric_tensor().inverse()

        H = [function('H1')(*x), function('H2')(*x), function('H3')(*x)]
        for i in range(1, 3+1):
            exc[i,0] = 0 # - H[i-1] / c # magnetic field
            # with the minus sign
        # D = [function('E1')(*x), function('E2')(*x), function('E3')(*x)]
        exc[1,2] = 0 # function('D12')(*x) # D[3-1] # dx∧dy
        # exc[2,3] = x[1]^2 * sin(x[2]) * g_inv[0, 0] * g_inv[1, 1] * function('epsilon')(*x) * function('Er')(*x)  # function('Er')(*x) # function('D23')(*x) # D[1-1] # dy∧dz
        ee = var("epsilon"); exc[2,3] = x[1]^2 * sin(x[2]) * g_inv[0, 0] * g_inv[1, 1] * ee * x[1] * function('Er')(*x)  
        exc[3,1] = 0 # function('Er')(*x) # function('D31')(*x) # D[2-1] # dz∧dx
        return exc

def compute_second_fundamental(spacetime, mode):
    x = spacetime.spherical_chart()
    U = spacetime.global_domain()
    surface = U.scalar_field(x[1]) # r = const.
    partials = spacetime.coordinate_frame()
    _dim = spacetime.dimension
    if mode == "orthonormal":
        frame = spacetime.compute_orthonormal_frame(partials=partials)
        normal_form = [frame[mu](surface) for mu in range(_dim)]
        eta = spacetime.sign_matrix
        C = spacetime.compute_structure_coefficients(frame)
        gamma = all_upper_coefficients(eta, C, g_inv = eta, domain = U)
    if mode == "coordinate": 
        g = spacetime.metric_tensor()
        gamma = g.christoffel_symbols()
        gamma = gamma[:]
        frame = partials
        normal_form = [frame[mu](surface) for mu in range(_dim)]
    K = second_fundamental_form(gamma, normal_form)
    return K

def print_second_fundamental_form(spacetime, mode, K):
    g = spacetime.metric_tensor()
    _dim = spacetime.dimension
    _ans = [[0 for i in range(_dim)] for j in range(_dim)]
    for i in range(_dim): 
        for j in range(_dim): 
            if K[i][j] != 0:
                print(i, j)
                print("second_fundamental_form K:")
                print(K[i][j].display())
                if mode == "orthonormal":
                    _ans[i][j] = K[i][j] / K[3][3]
                    print("K_[3][3] is set to be one", _ans[i][j].display())
                    
                if mode == "coordinate": 
                    _ans[i][j] = K[i][j] / g[i, j]
                    print("K[i][j] / g[i, j] = ", _ans[i][j].display())
                    print("metric g:")
                    print(g[i, j].display())
    print(g.display())
    return _ans
###########################################################
def set_const():
    _dim = 4
    M = Manifold(_dim, 'M', structure='Lorentzian')
    sm = SphericalMetric(M, f_r = "exp", g_r = "exp")
    return _dim, M, sm

def lhs(_dim, M, sm):
    g = sm.metric_tensor()
    R = g.ricci_scalar()
    Ric = g.ricci()
    print("g.display", g.display())
    print("g.ricci", Ric.display())
    print("g.ricci_scalar", R.display())

    # Set components in the coordinate frame e_X
    einstein_tensor = Ric.up(g, 1).copy() # M.tensor_field(1, 1, name='G', latex_name='G')
    for i, j in product(range(_dim), repeat=2): 
        einstein_tensor[i, j] -= kronecker_delta(i, j) * R / 2 # minus sign

    # Display the components (will show 1 on diagonals)
    einstein_tensor.display()
    for i, j in product(range(_dim), repeat=2): 
        if einstein_tensor[i, j] != 0:
            print(latex(einstein_tensor[i, j]))
    return einstein_tensor
def rhs(_dim, M, sm):
    F = sm.electromagnetic()
    H = sm.excitation()
    print("F.display()", F.display())
    print("H.display()", H.display())
    assume(0 < theta, theta < pi)
    # assume(sm.f_r.expr(), real)
    kT = sm.elemag_tensor11_FH(F, H)
    print("\n--- (kT)^0_0 の具体的な式 ---")
    # .expr() でシンボリックな式を取得
    for a, b in product(range(_dim), repeat=2):
        if kT[a, b] != 0:
            show(a, b)
            print(kT[a, b].expr())
            print(latex(kT[a, b].expr()))
    return kT
    
def main():
    _dim, M, sm = set_const()
    einstein_tensor = lhs(_dim, M, sm)
    kT = rhs(_dim, M, sm)
    return 0

if __name__ == '__main__':
    main()
    

###########################################################
