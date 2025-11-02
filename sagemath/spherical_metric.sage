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
        else:
            self.f_r = f_r
        if g_r is None:
            self.g_r = function('g')(r)
        else:
            self.g_r = g_r

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
# 使用例:
def main():
    M = Manifold(4, 'M', structure='Lorentzian')
    sm = SphericalMetric(M)
    g = sm.metric_tensor()
    g.display()
    return 0

if __name__ == '__main__':
    main()
###########################################################
