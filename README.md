# sac-cgt-convergence


editing


cgtbase.ipynb is main sym solver.

sactools.py and symcontools.py is library.

simulink/SAC_noise_ident_simu_3.slx is simulink SAC simulation file. it need initialization by plantset.m (Sm, Sp) and params/setting_3d.json.

CGT equation solver usage:
```python
pdim = 3; rdim = pdim # plant and refmodel dimention (same)
resultdict = ks.CGT_def_and_solve(pdim, rdim)
display(resultdict["kxkuEq"])
```


result:

```math
\begin{bmatrix} \\
a_{m3} b_{p2} - a_{m6} b_{p1} & a_{m3} b_{p3} - {a_{m}}_{9} {b_{p}}_{1} & {b_{p}}_{1} & 0 \\
{a_{m}}_{3} {b_{p}}_{3} - {a_{m}}_{9} {b_{p}}_{1} & {a_{m}}_{6} {b_{p}}_{3} - {a_{m}}_{9} {b_{p}}_{2} + {b_{p}}_{1} & {b_{p}}_{2} & 0 \\
{b_{p}}_{1} & {b_{p}}_{2} & {b_{p}}_{3} & 0 \\
{b_{m}}_{1} {b_{p}}_{2} - {b_{m}}_{2} {b_{p}}_{1} & {b_{m}}_{1} {b_{p}}_{3} - {b_{m}}_{3} {b_{p}}_{1} & 0 & {b_{p}}_{1} 
\end{bmatrix}

\begin{bmatrix}
k_{x11} \\
k_{x12} \\
k_{x13} \\
k_{u11}
\end{bmatrix}
=
\begin{bmatrix}
a_{m3} - a_{p3} \\
a_{m6} - a_{p6} \\
a_{m9} - a_{p9} \\
b_{m1}
\end{bmatrix}
```

