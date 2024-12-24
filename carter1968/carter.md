# metrics

Carter (1968) は時空の対称性を仮定して以下のような計量を研究した．

$ds^{2}=\dfrac{Z}{\Delta _{\lambda} }d\lambda ^{2}+\dfrac{Z}{\Delta _{\mu} }d\mu ^{2} + \dfrac{\Delta _{\mu} }{Z}\left[ P _{\lambda} d\psi -Q _{\lambda} d \chi \right] ^{2} -\dfrac{\Delta _{\lambda} }{Z}\left[ P_{\mu}d\psi -Q_{\mu}d\chi \right] ^{2}$

$Z=P_{\lambda }Q_{\mu}-P_{\mu }Q_{\lambda}$

Carter は電磁場を入れて目的である Hamilton--Jacobi 方程式の (変数) 分離可能性を議論した．

古典的 4元ポテンシャルは 1形式

$A = \dfrac{P_{\lambda }X_{\mu }+P_{\mu }X_{\lambda}}{Z}d\psi - \dfrac{Q_{\lambda }X_{\mu }+Q_{\mu }X_{\lambda}}{Z}d\chi $

のような形をしており，電磁場テンソルはその外微分

$F = dA$ (論文では × 2).

この電磁場テンソルと エネルギー・運動量テンソル の関係は
$4\pi T_{\mu \nu }= g^{\alpha \beta }F_{\mu \alpha }F_{\nu \beta } -\dfrac{1}{4}g_{\mu \nu }F_{\alpha \beta }F^{\alpha \beta }$

(ref. Gravitation p. 141, p. 471)

行列計算するには ($F_{\alpha \beta }F^{\alpha \beta } = F_{\alpha \beta  } g^{\alpha \alpha' } F_{\alpha' \beta' } g^{\beta' \beta } $ なので)

$4\pi (T_{\mu \nu })= (F_{\mu \alpha }) (g^{\alpha \beta }) (F_{\beta \nu })^{T} -\dfrac{1}{4}(g_{\mu \nu }) \cdot \operatorname{tr}[(F_{\alpha \beta })^{T} (g^{\alpha \alpha' }) (F_{\alpha' \beta' }) (g^{\beta' \beta })]$

あるいは

$4\pi T = ... $

を用いる．

## $A$

$ds^{2}=\left( \lambda ^{2}+\mu ^{2}\right) \left[ \dfrac{d\lambda ^{2}}{\Delta _{\lambda }}+\dfrac{d\mu ^{2}}{\Delta _{\mu} }\right] +\dfrac{\Delta_{\mu}\left[d \chi-\lambda^{2} d \psi\right]^{2}-\Delta_{\lambda}\left[d \chi+\mu^{2} d \psi\right]^{2}}{\lambda^{2}+\mu^{2}}$

* $Z=\lambda ^{2}+\mu ^{2}$
* $\Delta _{\lambda }=\dfrac{1}{3}\Lambda \lambda ^{4}+h\lambda ^{2}-2m\lambda +p+e^{2}$
* $\Delta _{\mu }=\dfrac{1}{3}\Lambda \mu ^{4}-h\mu ^{2}+2q\mu +p$
* $P\left( \lambda \right) =-\lambda ^{2}$
* $Q\left( \lambda \right) =-1$
* potential $A = e\left\{ \frac{\lambda \mu(\mu \cos \alpha+\lambda \sin \alpha)}{\lambda^{2}+\mu^{2}} d \psi+\frac{\lambda \cos \alpha-\mu \sin \alpha}{\lambda^{2}+\mu^{2}} d \chi\right\} $

Ricci


Field strength:

$F=dA$


## $\tilde{B}(+)$

$d s^{2}=\left(\lambda^{2}+l^{2}\right)\left\{\dfrac{d \lambda^{2}}{\Delta_{\lambda}}+\dfrac{d \mu^{2}}{\Delta_{\mu}}+\Delta_{\mu} d \psi^{2}\right\}-\dfrac{\Delta_{\lambda}[d \chi+2 l \mu d \psi]^{2}}{\lambda^{2}+l^{2}}$

* $Z=\lambda ^{2}+l ^{2}$
* $\Delta _{\lambda }=\Lambda \left( \dfrac{1}{3}\lambda ^{4}+2l^{2}\lambda ^{2}-l^{4}\right) + h\left( \lambda ^{2}-l^{2}\right) -2m\lambda +e^{2}$
* $\Delta _{\mu} =-h\mu ^{2}+2q\mu +p$
* $A=e\left\{\frac{\mu\left[2 l \lambda \cos \alpha+\left(\lambda^{2}+l^{2}\right) \sin \alpha\right]}{\lambda^{2}+l^{2}} d \psi+\frac{\lambda \cos \alpha}{\lambda^{2}+l^{2}} d \chi\right\}$


## $\tilde{B}(-)$

$d s^{2}=\left(\mu ^{2}+k^{2}\right)\left\{\dfrac{d \lambda^{2}}{\Delta_{\lambda}}+\dfrac{d \mu^{2}}{\Delta_{\mu}}-\Delta_{\lambda} d \psi^{2}\right\}+\dfrac{\Delta_{\lambda}[d \chi-2 k \lambda d \psi]^{2}}{\mu ^{2}+k^{2}}$

* $Z=\mu ^{2}+k^{2}$
* $\Delta _{\lambda} =h\lambda ^{2}-2m\lambda +n$
* $\Delta _{\mu }=\Lambda \left( \dfrac{1}{3}\mu ^{4}+2k^{2}\mu ^{2}-k^{4}\right) - h\left( \mu ^{2}-k^{2}\right) +2q\mu -e^{2}$
* $A=e\left\{\frac{\lambda\left[\left(\mu^{2}+k^{2}\right) \cos \alpha + 2 k \mu \sin \alpha\right]}{\mu^{2}+k^{2}} d \psi+\frac{\mu \sin \alpha}{\mu^{2}+k^{2}} d \chi\right\}$

## $D$

$ds^{2}=\dfrac{d\lambda^{2}}{\Delta _{\lambda}}+\dfrac{d\mu ^{2}}{\Delta _{\mu} }+\Delta _{\mu }d\chi^{2}-\Delta _{\lambda}d\psi ^{2}$

* $Z=1$
* $\Delta_{\lambda}=\left(\Lambda+e^{2}\right) \lambda^{2}-2 m \lambda+n $
* $\Delta_{\mu}=\left(\Lambda-e^{2}\right) \mu^{2}+2 q \mu+p$
* $P\left( \lambda \right) =0$
* $Q\left( \lambda \right) =-1$
* $A=e\{\lambda \cos \alpha d \psi-\mu \sin \alpha d \chi\}$

電磁場テンソルの計算

$F/2=dA=e\{\cos \alpha d\lambda \wedge d\psi -\sin \alpha d\mu \wedge d\chi\}$

だから $F_{\lambda \psi }=e\cos \alpha, F_{\mu \chi}=-e\sin \alpha  $

$\left( \lambda ,\mu ,\chi,\psi\right) $ の順番で書くと


## $\tilde{A}$

$ds^{2}=\left[ \left( c\lambda \cos \gamma +k\right) ^{2}+ \left( c {\mu }\sin \gamma +l\right) ^{2}\right] \left\{ \dfrac{d\lambda ^{2}}{\Delta _{\lambda }}+\dfrac{d\mu ^{2}}{\Delta _{\mu }}\right\} $

$+\dfrac{\Delta \left\{ \sin \gamma d\chi-[\left( c^{2}\lambda ^{2}+k^{2}+\lambda ^{2}\right) \cos \gamma + 2ck\lambda ] d\psi \right\} ^{2}}{\left( c\lambda \cos \gamma +k\right) ^{2} +\left( c {\mu }\sin \gamma +l\right) ^{2}}$
$-\dfrac{\Delta \left\{ \cos \gamma d\chi+[\left( c^{2}\mu ^{2}+k^{2}+\lambda ^{2}\right) \sin \gamma + 2cl\mu ] d\psi \right\} ^{2}}{\left( c\lambda \cos \gamma +k\right) ^{2} +\left( c {\mu }\sin \gamma +l\right) ^{2}}$

where 

$\Delta _{\lambda }=\dfrac{1}{3}\Lambda c^{2}\lambda ^{4}\left( \cos \gamma \right) ^{2} +\dfrac{4}{3}\Lambda ck\lambda ^{3}\cos \gamma +\left( 2\Lambda l^{2}+h\right) \lambda ^{2} -2m\lambda +n$

$\Delta _{\mu }=\dfrac{1}{3}\Lambda c^{2}\mu ^{4}\left( \sin \gamma \right) ^{2} +\dfrac{4}{3}\Lambda cl\mu ^{3}\sin \gamma +\left( 2\Lambda k^{2}-h\right) \mu ^{2} +2q\mu +p$

$e^{2} = \Lambda \left( k^{4}-l^{4}\right) +h\left( k^{2}+l^{2}\right) + 2c\left( km\cos \gamma +lq\sin \gamma \right) +c^{2}\left( n\left( \cos \gamma \right) ^{2}-p\left( \sin \gamma \right) ^{2}\right) $