# Details of the Decay Amplitude


[1. The angular distribution $W_l$](#theangulardistribution)  
[2. The barrier factor $B_l$](#thebarrierfactor)  
[3. The propagators](#thepropagators)  
  [3.1 Breit-Wigner](#breitwigner)  
  [3.2 Gounaris-Sakurai](#gounarissakurai)  
  [3.3 LASS](#lass)  
  [3.4 Generalized LASS](#generalizedlass)  
  [3.5 K-Matrix](#kmatrix)  
  [3.6 Flatté](#flatte)  


### 1. The angular distribution $W_l$ <a name="theangulardistribution"></a> 

The angular distribution $ W_l(m^2_{ab},m^2_{bc}) $ follows the Zemach formalism[^Zemach:1963bc][^Zemach:1965ycj][^Filippini:1995yc]. For the $ P\to cr $, $ r\to ab $ case, $ W_l(m^2_{ab},m^2_{bc}) = Z_l(m^2_{ab},m^2_{bc})$, and: 


$$
\begin{align}
Z_0 & = 1 \quad ,\nonumber\\  
Z_1 & = m^2_{bc}-m^2_{ac}+\frac{(m_P-m^2_{c})(m^2_{a}-m^2_{b})}{m^2_{ab}} \quad ,\nonumber\\ 
Z_2 & = \Big[m^2_{bc}-m^2_{ac}+\frac{(m_P-m^2_{c})(m^2_{a}-m^2_{b})}{m^2_{ab}}\Big]^2 \nonumber \\
    & - \frac{1}{3}\Big[m^2_{ab}-2(m^2_{P}+m^2_{c})+\frac{m^2_{P}-m^2_{c}}{m^2_{ab}}\Big]\Big[m^2_{ab}-2(m^2_{a}+m^2_{b})+\frac{m^2_{a}-m^2_{b}}{m^2_{ab}}\Big]  \quad . \nonumber
\end{align}
$$

For the decay amplitude component $ P\to ar $, $ r\to bc $,  $ W_l(m^2_{ab},m^2_{bc}) = Z_l(m^2_{bc},m^2_{ac})$. Similarly, one obtains the $ W_l(m^2_{ab},m^2_{bc}) $ for $ P\to ar $, $ r\to bc $.  The $ m^2_{ac} $ is obtained by the four momentum conservation. And there is a "helicity" switch in DAFNE for the helicity formalism of angular distribution. After turning on the "helicity" option, the $ m^2_{ab} $ in denominator will change into the resonance mass squared $ m^2_r $. 



### 2. The barrier factor $B_l$ <a name="thebarrierfactor"></a> 
The barrier factors $ B_l^{Pc}(|\vec{p}|,|\vec{p}_r|) $ and $ B_l^{ab}(|\vec{q}|,|\vec{q}_r|) $, are related to the momentum of $ c $ in the mother particle rest frame ($ \vec{p} $) and the momentum of $ a $ in the $ ab $ system rest frame ($ \vec{q} $). The $ \vec{p}_r $ ($ \vec{q}_r $) represents the $ \vec{p} $ ($ \vec{q} $) when $ m_{ab} $ is right at the resonance mass $ m_r $. The form $ B_l(|\vec{k}|,|\vec{k}_r|) $ are the Blatt-Weisskopf centrifugal barrier factors[^Aston:1987ir][^VonHippel:1972fg][^Blatt:1952ije]:
$$
\begin{align}
B_0 &= 1 \quad , \nonumber \\
B_1 &= \frac{\sqrt{z_r^2+1}}{\sqrt{z^2+1}} \quad , \nonumber \\
B_2 &= \frac{\sqrt{z_r^4+3z_r^2+9}}{\sqrt{z^4+3z^2+9}} \quad ,
\end{align}
$$
where $ z = R|\vec{k}|$ and $ z_r = R|\vec{k}_r| $, with Blatt-Weisskopf effective radius $ R $. There are multiple definition of barrier factors, and this form is used by ref[^BaBar:2010nhz][^cfit][^Jordi:thesis] and the Hydra framework[^hydra]. 



### 3. The propagators <a name="thepropagators"></a> 
#### 3.1 Breit-Wigner <a name="breitwigner"></a> 

The Breit-Wigner function is an approximation for the shape of narrow resonance. In the code, we follow the relativistic Breit-Wigner: 
$$
\begin{equation}
T_{\rm BW}(m_{ab}) = \frac{1}{m^2_r-m^2_{ab}-im_r\Gamma_r(m_{ab})} \quad ,
\end{equation}
$$
where the running width
$$
\begin{equation}
\Gamma_r(m_{ab})=\Gamma^0_r\bigg(\frac{|\vec{q}|}{|\vec{q_r}|}\bigg)^{2l+1}\bigg(\frac{m_r}{m_{ab}}\bigg)\big(B_l^{ab}(|\vec{q}|,|\vec{q}_r|)\big)^2 \quad .
\end{equation}
$$
#### 3.2 Gounaris-Sakurai <a name="gounarissakurai"></a> 
For the $ \pi\pi $ P-Wave, the Gounaris-Sakurai function[^Gounaris:1968mw] describes better then the Breit-Wigner function. Mathematically, Gounaris-Sakurai function has better analyticity. We follow the form: 
$$
\begin{equation}
T_{\rm GS}(m_{ab}) = \frac{1+\Gamma^0_rd/m_r}{m^2_r-m^2_{ab}+f(m^2_{ab})-im_r\Gamma_r(m_{ab})} \quad ,
\end{equation}
$$
where
$$
\begin{equation}
d = \frac{3(m^2_a+m^2_b)}{2\pi|\vec{q}_r|^2}\ln(\frac{m_r+2|\vec{q}_r|}{m_a+m_b}) + \frac{m_r}{2\pi|\vec{q}_r|} - \frac{(m^2_a+m^2_b)m_r}{2\pi |\vec{q}_r|^3}
\end{equation}
$$
$$
\begin{equation}
f = \frac{\Gamma_r^0 m_r^2}{|\vec{q}_r|}\Big\{\frac{|\vec{q}|^2}{|\vec{q}_r|^2}[h(m^2_{ab})-h(m^2_r)]+(m^2_r-m^2_{ab})h^\prime(m^2_r)\Big\} \quad ,
\end{equation}
$$
and
$$
\begin{align}
h(m^2) & = \frac{2|\vec{q}|}{\pi m}\ln\bigg(\frac{m_{ab}+2|\vec{q}|}{m_a + m_b}\bigg) \\ 
h^\prime(m_r^2) & = \frac{dh}{d(m^2)}\bigg|_{m^2=m^2_r} = \bigg(\frac{1}{8|\vec{q}_r|^2} - \frac{1}{2m^2_r}\bigg)h(m^2_r) + 	\frac{1}{2\pi m^2_r} \quad .
\end{align}
$$

#### 3.3 LASS <a name="lass"></a> 
For the $ K\pi $ S-Wave resonance, $ K^*(1430) $, LASS model describes the shape better then Breit-Wigner. This model is a combination of a resonance part and a non-resonance part. We follow the form[^Lupton:2016iaz]:
$$
\begin{equation}
T_{\rm LASS}(m_{ab}) = f_r\Big(\frac{m_{ab}}{m_r}\Big)sin(\delta_R+\delta_B)e^{i(\delta_R+\delta_B)} \quad .
\end{equation}
$$
The empirical form factor
$$
\begin{equation}
	f_r(x) = \exp(b_1x + b_2x^2 + b_3x^3)
\end{equation}
$$
The phases
$$
\begin{align}
\delta_R & = \tan^{-1}\bigg(\frac{m_r\Gamma_r^0}{m_r^2-m^2_{ab}}\bigg) \nonumber \quad ,\\
\delta_B & = \cot^{-1}\bigg(\frac{1}{a|\vec{q}|}+\frac{r|\vec{q}|}{2}\bigg)\nonumber \quad .
\end{align}
$$

#### 3.4 Generalized LASS <a name="generalizedlass"></a> 
With more parameters, Generalized LASS can describe the $ K^*(1430) $ from decays better then LASS. We follow the form: 
$$
\begin{align}
T_{\rm GLASS} = & Re^{i\phi_R+2i\phi_B}\Big(\frac{q \cot\delta_B+iq}{q \cot\delta_B-iq}\Big)\frac{m_r}{2|\vec{q}_r|}\Big(\frac{m_r\Gamma_r^0}{m^2_r - m^2_{ab} -im_r\Gamma(m_{ab})}\Big) \nonumber \\
& + B\frac{m_{ab}}{2}e^{i\phi_B}\frac{\cos\phi_B+\sin\phi_B\cot\delta_B}{q\cot\delta_B-i|\vec{q}|}
\end{align}
$$
The Generalized LASS reduces to LASS when setting $ B=R=1 $ and $ \phi_B = \phi_R = 0 $.

#### 3.5 K-Matrix <a name="kmatrix"></a> 
In the case where multiple channels coupling becomes important, the K-Matrix is more effective than the sum of isolated resonances, for example, when describing the $ \pi\pi $ S-Wave. The K-Matrix $ \hat{K} $ parameterizes the scattering processes with multiple resonances that indicated with poles. The scattering process is usually parameterized as
$$
\begin{equation}
(1-i\hat{K}\rho)^{-1}\hat{K} = \hat{K}-\hat{K}^2\rho+\hat{K}^3\rho^2+\cdots
\end{equation}
$$
This form can be generalized to describe resonances in meson decays:
$$
\begin{equation}
T_{\rm K-Matrix}=(1-i\hat{K}\rho)^{-1}\hat{P} \quad .
\end{equation}
$$
The production vector $ \hat{P} $ parameterizes the resonances production in meson decay, while the $ (1-i\hat{K}\rho)^{-1} $ part describe the decays of the resonances. For $ \hat{K} $, we follow the form 
$$
\begin{align}
\hat{K}_{ij} = \bigg(\sum_{\alpha}\frac{g^0_{\alpha i}g^0_{\alpha j}}{m^2_{\alpha}-m^2_{ab}}+f^{sc}_{ij}\frac{1{\rm GeV^2}-s^{sc}_0}{m^2_{ab}-s^{sc}_0}\bigg)\bigg[\frac{1{\rm GeV^2}-s^A_0}{m^2_{ab}-s^A_0}\Big(m^2_{ab}-\frac{s_A m_a m_b}{2}\Big)\bigg] \quad .
\end{align}
$$
$ g_{\alpha i} $ is the real constant describing the coupling between the resonance $ \alpha $ and the channel $ i $. In our code, $ i $ loops over $ \pi^+\pi^- $, $ K\bar{K} $, $ 4\pi $, $ \eta\eta $ and $ \eta\eta^{\prime} $. $ f^{sc}_{ij} $ and $ s^{sc}_0 $ describe a smooth part that do not related to the resonance $ \alpha $. The square brackets envelop the Adler zero term[^Adler:1964um]. The phase space factor, with analytical continuation, follows the definition in ref[^Anisovich:2002ij]. The $ \hat{P} $ has a form similar to $ \hat{K} $ 
$$
\begin{equation}
\hat{P}_j = \sum_{\alpha}\frac{\beta^0_\alpha g^0_{\alpha j}}{m^2_{\alpha}-m^2_{ab}}+f^{pr}_{kj}\frac{1{\rm GeV^2}-s^{pr}_0}{m^2_{ab}-s^{pr}_0} \quad .
\end{equation}
$$
$ \beta^0_\alpha $ is the complex coefficient for the resonance $ \alpha $. The complex parameter $ f^{pr}_{kj} $ and the real $ s^{pr}_0 $ describe the smooth part, for the reconstructed channel $ k $. Our code sets $ k = \pi\pi $ as default. The $ g_{\alpha i} $, $ m_\alpha $, $ f^{sc}_{ij} $, $ s^{sc}_0 $, $ s^{A}_0 $ and $ s_A $ are fixed to the values from the scattering experiment[^Jordi:thesis][^Anisovich:2002ij]. Currently, only the $ i = \pi\pi $ and $ j = \pi\pi $ parts of $ f^{sc}_{ij} $ is filled, while other parts are left zero.

#### 3.6 Flatté <a name="flatte"></a> 
When there is only one resonance and two coupling channels, the K-Matrix reduced to the Flatté model[^Ayik:1976mqe]. We follow the form:
$$
\begin{equation}
T_\text{Flatté} = \frac{m_r\Gamma^0_r\gamma_{r1}^2}{m^2_r-m^2_{ab}-im_r\Gamma^0_r(\gamma_{r1}^2\rho_1+\gamma_{r2}^2\rho_2)}
\end{equation}
$$
The $ \gamma_{r1} $ and $ \gamma_{r2} $ are the fit parameters related to the coupling between the channel 1,2 and the resonance. $ \rho_1 $ and $ \rho_2 $ are the phase space factor of the channel 1 and 2. Channel 1 is the reconstructed channel, and has definite masses for the two final state particles. The masses of that in channel 2 are included in fit parameters. 



[^Zemach:1963bc]:  C.~Zemach,  Three pion decays of unstable particles, Phys. Rev. **133** (1964), B1201. doi:10.1103/PhysRev.133.B1201
[^Zemach:1965ycj]:  C.~Zemach, Use of angular momentum tensors, Phys. Rev. **140** (1965), B97-B108. doi:10.1103/PhysRev.140.B97
[^Filippini:1995yc]:  V.~Filippini, A.~Fontana and A.~Rotondi, Covariant spin tensors in meson spectroscopy, Phys. Rev. D **51** (1995), 2247-2261. doi:10.1103/PhysRevD.51.2247
[^Aston:1987ir]:  D.~Aston, N.~Awaji, T.~Bienz, F.~Bird, J.~D'Amore, W.~Dunwoodie, R.~Endorf, K.~Fujii, H.~Hayashi and S.~Iwata, *et al.*, A Study of $K^- \pi^+$ scattering in the reaction $K^- p \to K^- \pi^+ n$ at $11$ ${\rm GeV}/c$, Nucl. Phys. B **296** (1988), 493-526. doi:10.1016/0550-3213(88)90028-4
[^VonHippel:1972fg]:  F.~Von Hippel and C.~Quigg, Centrifugal-barrier effects in resonance partial decay widths, shapes, and production amplitudes, Phys. Rev. D **5** (1972), 624-638. doi:10.1103/PhysRevD.5.624
[^Blatt:1952ije]:  J.~M.~Blatt and V.~F.~Weisskopf, Theoretical nuclear physics, Springer, 1952, ISBN 978-0-471-08019-0. doi:10.1007/978-1-4612-9959-2
[^BaBar:2010nhz]:  P.~del Amo Sanchez *et al.*, Measurement of $D^0$-$\overline{D}^0$ mixing parameters using $D0 \to K_S^0 \pi^+ \pi^-$ and $D^0 \to K_S^0 K^+ K^-$ decays, Phys. Rev. Lett. **105** (2010), 081803. arXiv:1004.5053, doi:10.1103/PhysRevLett.105.081803
[^cfit]:  J.~G.~Tico, cfit. https://github.com/cfit/cfit
[^Jordi:thesis]:  J.~G.~Tico, Measurement of the neutral D meson mixing parameters at the BaBar experiment.
[^hydra]:  Alves Junior, A A, MultithreadCorner/Hydra. doi:10.5281/zenodo.1206261
[^Gounaris:1968mw]:  G.~J.~Gounaris and J.~J.~Sakurai, Finite width corrections to the vector meson dominance prediction for $\rho \to e^+ e^-$, Phys. Rev. Lett. \textbf{21} (1968), 244-247. doi:10.1103/PhysRevLett.21.244
[^Lupton:2016iaz]:  O.~Lupton, Studies of $\mathrm{D}^0\to\mathrm{K}^0_{\scriptscriptstyle\rm S}\mathrm{h}^{+}\mathrm{h}'^{-}$ decays at the LHCb experiment, CERN-THESIS-2016-156.
[^Adler:1964um]:  S.~L.~Adler, Consistency conditions on the strong interactions implied by a partially conserved axial vector current, Phys. Rev. **137** (1965), B1022-B1033. doi:10.1103/PhysRev.137.B1022
[^Anisovich:2002ij]:  V.~V.~Anisovich and A.~V.~Sarantsev, K matrix analysis of the ($I J^{PC} = 00^{++}$)-wave in the mass region below 1900 MeV, Eur. Phys. J. A **16** (2003), 229-258. arXiv:hep-ph/0204328, doi:10.1140/epja/i2002-10068-x
[^Ayik:1976mqe]:  S.~Ayik, Generalized master equations for heavy-ion collisions, Phys. Lett. B **63** (1976), 22-24. doi:10.1016/0370-2693(76)90458-5