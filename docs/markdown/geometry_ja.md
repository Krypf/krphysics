# 目的
<!-- setup -->
$M$: $m$ 次元可微分多様体  
$N$: $n$ 次元可微分多様体  

<!-- end setup -->

外在曲率を2次共変テンソルとして定義する．

# 準備

## 埋め込み

<!-- def. -->
**定義.** 位相空間 $\left( X,T\right)$ から位相空間 $\left( Y,S\right)$ への**埋め込み**とは，**$X$ から $Y$ の部分空間への同相写像**である．つまり，$X$ から $Y$ への埋め込みとは，写像 $f: X \rightarrow Y$ であって $(f\left[ X \right], S_{f\left[ X \right]})$ が $\left( X,T\right)$ と同相となるものである．
<!-- end def. -->
<!-- https://en.wikipedia.org/wiki/Embedding -->
<!-- https://www.rimath.saitama-u.ac.jp/lab.jp/Fukui/lectures/Set_Topsp.pdf -->
特に，以下の性質を満たす．

1. 埋め込みは連続
2. 埋め込みは単射

<!-- ex. -->
**例.** 
$X=\lbrace 1,2\rbrace , Y = \lbrace 1,2,3\rbrace$ および $T=2^{X} = \lbrace \lbrace \rbrace , \lbrace 1\rbrace, \lbrace 1,2\rbrace \rbrace$, $S=\lbrace \lbrace \rbrace , \lbrace 1\rbrace, \lbrace 1,2\rbrace ,\lbrace 3,1 \rbrace, Y \rbrace$ とする．

$(f(1), f(2)) = (1, 2)$ で定まる写像 $f$ は $(X,T)$ から $(Y,S)$ への埋め込みである．

一方，$(f(1), f(2)) = (2, 1)$ で定まる異なる写像 $f$ は，$f^{-1}[\{1\}] = \{2\}$ を満たし，連続写像ではないので埋め込みでもない．
<!-- end ex. -->

<!-- ex. -->
**例.** $n$ 次元 実数 空間 $\mathbb{R}^{n}$ ($n \in \mathbb{N} $) は $(n + 1)$ 次元 実数 空間 $\mathbb{R}^{n+1}$ において，超平面 $P = \lbrace \left( x_{j}\right) _{j=1}^{n+1}| x_{n+1}=0 \text{ and } ∀ j \in \{1,...,n\},  x_{j}\in \mathbb{R} \rbrace$ として埋め込まれる．ただし，ここで埋め込みは 包含写像 $\iota\colon P \to \mathbb{R}^{n+1} $.
<!-- end ex. -->

<!-- ex. -->
**例.** 直前の例で 次元の組 $(n, n+1)$ を 一般に $(m,n) \, (m \le n )$ としても同様に埋め込みが成り立つ．
<!-- end ex. -->

## はめ込み

<!-- setup -->
$M$: $m$ 次元 $C^{r}$ 級多様体  
$N$: $n$ 次元 $C^{r}$ 級多様体  
$f\colon M \to N$: $C^{r}$ 級写像
<!-- end setup -->

<!-- def. -->
**定義.** $M$ 上で $f $ の微分が単射であるとき，写像 $f$ を $M$ から $N$ への **はめ込み** という． すなわち，はめ込み $f$ について，

$\forall p \in M, \, (df)(p) \equiv (df)_{p} \colon T_{p}M\rightarrow T_{f\left( p\right) }N$ は単射である．

<!-- end def. -->

<!-- def. -->
**定義.** はめ込みであり，かつ位相空間の間の (微分同相な) 埋め込みでもある写像を，(可微分多様体の間の) **埋め込み**と呼ぶ．
<!-- end def. -->

ただし，埋め込みを単に単射なはめ込みと定義することもある．
<!-- http://www2.itc.kansai-u.ac.jp/~afujioka/2012-2016/2014/g4/141113g4.pdf -->
<!-- Nakahara -->

<!-- def. -->
**定義.** 包含写像 $\iota \colon M \to N$ が埋め込みであるとき，$M$ は $N$ の **部分多様体** であるという．
<!-- end def. -->

<!-- ex. -->
**例.** $(M,N) = (\mathbb{R}^{m}, \mathbb{R}^{n})$ (ただし $M \subsetneq N$ とみなす) について，包含写像 $\iota \colon M \to N$ は 接空間の間の包含写像 $\partial _{1} \mapsto \partial _{1}, \dots, \partial _{m} \mapsto \partial _{m}$ を定め，これは単射なので $M$ は $N$ の部分多様体といえる．
<!-- end ex. -->

# 法バンドル

<!-- setup -->
* $\left( N,g_{N}\right)$:  $n$ 次元 一般 Riemann 多様体 (とりわけ Lorentz 多様体)  
* $M$ を $N$ の $m$ 次元部分多様体とする．
* $M$ の アトラス（座標近傍系）を $\mathcal{A}$, $N$ のアトラスを $\mathcal{B}$ とする．
<!-- end setup -->

## 誘導計量

$N$ の計量の記号を簡単に $g \equiv g_{N}$ と省略する．${\left( x,U\right)} \in \mathcal{A}, {\left( y,V\right)} \in \mathcal{B}$ とすると $M$ 上の $(N, g)$ からの誘導計量は

$g_{M}=\dfrac{\partial y^{k}}{\partial x^{i}}\dfrac{\partial y^{l}}{\partial x^{j}} g_{kl} dx^{i}\otimes dx^{j}$ (引数の点 $p \in U \cap V$ を省略した)
<!-- arai p. 292 -->

以降 $g_{M} \equiv h$ とする．すると 

$h_{ij} (p) = \dfrac{\partial y^{k}}{\partial x^{i}}(p) \dfrac{\partial y^{l}}{\partial x^{j}}(p)  g_{kl}(p)$.

## 接空間

次に $N$ 上の接空間 $(T_{p}N, g(p))$ を $M$ 上の接空間 $(T_{p}M, g_{M}(p))$ として，$T_{p}N$ の $T_{p}M$ に関する直交補空間を考える．

1. $e_{a}\left( p\right) =e_{a}^{j}\left( p\right) \dfrac{\partial }{\partial x^{j}}\left( p\right) \in T_{p}M$ 

2. $n_{b}(p) = n_{b} ^{j}(p) \dfrac{\partial }{\partial y^{j}} \left( p\right) \in T_{p} N  \, (j\in \lbrace 1,\ldots ,n\rbrace,\, b\in \lbrace 1,\ldots ,n-m\rbrace  ) $

以下，点の引数を省略する．$\dfrac{\partial }{\partial x^{j}}=\dfrac{\partial y^{k}}{\partial x^{j}}\dfrac{\partial }{\partial y^{k}}$ より $e_{a}=e_{a}^{j}\dfrac{\partial }{\partial x^{j}}=e_{a}^{j}\dfrac{\partial y^{k}}{\partial x^{j}}\dfrac{\partial }{\partial y^{k}}$ なので 直交条件 $g_{N}\left( e_{a},n_{b}\right) =0$ より，

$g_{\mu \nu }dy^{\mu }\otimes dy^{\nu }\left( e_{a}^{i}\dfrac{\partial y^{k}}{\partial x^{i}}\dfrac{\partial }{\partial y^{k}},n_{b}^{j}\dfrac{\partial }{\partial y^{j}}\right) = g_{\mu \nu }\left( p\right) e_{a}^{i}\left( p\right) \dfrac{\partial y^{\mu}}{\partial x^{i}}\left( p\right) n_{b}^{\nu}\left( p\right) = 0$

規格化条件より $g_{\mu \nu }\left( p\right) n_{b}^{\mu }\left( p\right) n_{b}^{\nu }\left( p\right) = \pm 1 $?

# 超曲面

* $\Sigma$: $N$ 上の $C^{r}$ 級超曲面  
* $D \coloneqq (n - 1)$ とおく．  

* $\Sigma$ は $D$ 次元部分多様体なので 座標近傍系 $\mathcal{A}$ が取れる．

