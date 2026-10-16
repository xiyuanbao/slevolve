# 一阶（degree-1）负荷 Love 数：参考系与配置法初始化缺陷

## 1. 一阶参考系
在一阶，地球的刚体平移是一个自由模态，因此必须先选定原点，h_1、l_1、k_1 才有确定的值。改变原点时，h_1、l_1、k_1 会同时加上同一个常数。常用约定有三种：
- **CM**（地球 + 负荷的质心）：k_1 = -1；
- **CE**（固体地球的质心）：k_1 = 0；
- **CF**（地表形状中心）：h_1 + 2 l_1 = 0。

`lovenum/` 中的所有 Love 数文件的弹性与流体极限都在 CE 系（k_1 = 0）（Case 1 文件的一阶采样值见第 2 节），SLEvolve 不做任何参考系转换：一阶位移与大地水准面输出处于 Love 数文件的参考系中。h_1 - k_1、l_1 - k_1 等差值与参考系无关；海平面取决于 1 + k - h，因此与参考系无关。

## 2. 一阶配置法初始化缺陷
配置法（Peltier 1985；Mitrovica and Peltier 1992）在固定极点 s_i 上拟合 sum_j a_j/(s_i + s_j) = X(s_i) - X_el，其中 X_el 是单独计算的弹性极限（s -> 无穷）。按构造，左边只表示松弛部分：s -> 无穷时为零。因此该拟合假设采样值 X(s_i) 与极限 X_el 处于同一一阶参考系。

在以往的库文件中（Austermann and Mitrovica 2015 的 Love 数库，Case 2 文件以前取自该库），一阶采样值在 CM 系中初始化，而极限值在 CE 系中。通俗地说：每个采样值都带有 h、l、k 中相同的常数偏差 -1，而衰减项 a_j/(s + s_j) 之和无法表示常数，于是拟合用很大的正负交替留数来补偿；一阶松弛曲线 h_1(t)、l_1(t)、k_1(t) 因此偏差同一个随时间变化的量。`lovenum/` 中的 Case 2 与 Case 3 文件在计算时采样值与极限值处于同一参考系，消除了该缺陷（其一阶 k 留数严格为 0）。Case 1 文件（`prem_512.l60K2C.sum18p6.dum19p2.tz19p4.lm22`，随 Han et al. 2022 分发，原样收录）带有同样的一阶缺陷，对海平面与地形没有影响。

由于同一误差同时进入 h、l、k，它在 h - k 中抵消：海平面与地形不受影响（一阶 h - k（无量纲 Love 数）的变化不超过 3e-5）。单独输出的一阶场（径向位移、大地水准面与水平位移）有差别。

## 参考文献
- Austermann, J. and Mitrovica, J. X.: Calculating gravitationally self-consistent sea level changes driven by dynamic topography, Geophys. J. Int., 203, 1909-1922, 2015.
- Han, H. K., Gomez, N., and Wan, J. X. W.: Geosci. Model Dev., 15, 1355-1373, 2022; data: doi:10.17605/OSF.IO/8PTFM.
- Mitrovica, J. X. and Peltier, W. R.: A comparison of methods for the inversion of viscoelastic relaxation spectra, Geophys. J. Int., 108, 410-414, 1992.
- Peltier, W. R.: The LAGEOS constraint on deep mantle viscosity: Results from a new normal mode method for the inversion of viscoelastic relaxation spectra, J. Geophys. Res., 90, 9411-9421, 1985.
