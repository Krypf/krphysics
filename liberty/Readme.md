# Topological Poset Betti Number Calculator

有限のアレクサンドロフ位相空間（Alexandrov Topology）および半順序集合（Poset）から順序複体（Order Complex）を構成し、その空間の位相不変量である **ベッチ数（Betti numbers, $b_k$）** を正確かつ高速に計算する Julia パッケージです。

## 1. 理論的背景

### アレクサンドロフ位相と半順序集合 (Poset)
有限位相空間では、開集合系と要素間の前順序関係（Preorder）が完全に一対一対応します。本プログラムでは、事象 $x, y$ の間に以下の特殊化前順序（Specialization Preorder）を定義します：
$$x \le y \iff \forall U \in \mathcal{O}, (y \in U \implies x \in U)$$

### 順序複体 (Order Complex)
順序集合の構造を幾何学的な単体複体（Simplicial Complex）に変換します。
*   **0-単体（頂点）**: 集合の各要素。
*   **$k$-単体（面）**: 長さ $k$ の全順序の鎖（Chains） $x_0 < x_1 < \dots < x_k.

### ベッチ数 (Betti Numbers)
ベッチ数は空間のトポロジー的な「穴」の数を表す不変量です。
*   **$b_0$**: 連結成分の数（独立した島の数）。
*   **$b_1$**: 1次元の穴（ループ・円環状の空洞）の数。
*   **$b_2$**: 2次元の穴（立体的な空洞）の数。

代数的には、$k$ 次ホモロジー群 $H_k = \ker \delta_k / \text{im } \delta_{k+1}$ のランクとして計算されます。
$$b_k = \vert{}C_k\vert{} - \text{rank}(\partial_k) - \text{rank}(\partial_{k+1})$$
ここで、$\partial_k$ は境界写像（Boundary Operator）を行列化したものです。

## 2. 特徴・最適化手法

1.  **非 $T_0$ 空間（Kolmogorov 商空間）の自動処理**:
    推移的な密着関係（$x \le y$ かつ $y \le x$）を自動検出し、強連結成分を1つの同値類として縮約（商空間化）することで、境界計算時の無限ループを防止します。
2.  **LinearAlgebra　による計算**:
    Juliaの標準ライブラリ `LinearAlgebra` を使用し、浮動小数点行列としてのランク計算を高速に行います。

## 3. 動作要件

*   **Julia** (v1.6以上を推奨)
*   標準パッケージ: `LinearAlgebra`, `Printf`

## 4. 使い方

モジュールを読み込み、全体集合（Universe）と開集合系のリスト（Setの配列）を渡すだけで計算が完了します。

```julia
using .PosetBetti

# 定義: 空間 X = {0, 1, 2, 3}
universe = [0, 1, 2, 3]

# 例: クラウン構造 (分散型社会モデル)
# オープンセット (自由の器) の定義
open_sets = [
    Set{Int}(), Set([0]), Set([1]), Set([0, 1]), 
    Set([2]), Set([0, 2]), Set([0, 1, 2])
]

# インスタンスの作成と計算
calc = PosetBettiCalculator(universe, open_sets)
betti_numbers = compute_all_betti_numbers(calc)

println(betti_numbers)
# 期待される出力: Dict(0 => 1, 1 => 1, 2 => 0)
# b_1 = 1 により、権力の空白(穴)が存在する分散調和構造であることがわかります。