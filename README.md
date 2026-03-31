# sac-cgt-convergence


editing


cgtbase.ipynb is main sym solver.

knw_symtools.py and symcontools.py is library.

SAC_noise_ident_simu_3.slx is simulink SAC simulation file. it need initialization by plantset.m (Sm, Sp) and params/setting_3d.json.

CGT equation solver usage:
```python
pdim = 3; rdim = pdim # plant and refmodel dimention (same)
resultdict = ks.CGT_def_and_solve(pdim, rdim)
display(resultdict["kxkuEq"])
```

ddd

$$
\begin{array}
{a_{m}}_{3} {b_{p}}_{2} - {a_{m}}_{6} {b_{p}}_{1} & {a_{m}}_{3} {b_{p}}_{3} - {a_{m}}_{9} {b_{p}}_{1} & {b_{p}}_{1} & 0 \\
{a_{m}}_{3} {b_{p}}_{3} - {a_{m}}_{9} {b_{p}}_{1} & {a_{m}}_{6} {b_{p}}_{3} - {a_{m}}_{9} {b_{p}}_{2} + {b_{p}}_{1} & {b_{p}}_{2} & 0 \\
{b_{p}}_{1} & {b_{p}}_{2} & {b_{p}}_{3} & 0 \\
{b_{m}}_{1} {b_{p}}_{2} - {b_{m}}_{2} {b_{p}}_{1} & {b_{m}}_{1} {b_{p}}_{3} - {b_{m}}_{3} {b_{p}}_{1} & 0 & {b_{p}}_{1}
\end{array} \quad \begin{array} {k_x}_{11} \\
{k_x}_{12} \\
{k_x}_{13} \\
{k_u}_{11} \end{array}
$$

result:  

```math
\left(
\left[\begin{matrix}
{a_{m}}_{3} {b_{p}}_{2} - {a_{m}}_{6} {b_{p}}_{1} & {a_{m}}_{3} {b_{p}}_{3} - {a_{m}}_{9} {b_{p}}_{1} & {b_{p}}_{1} & 0 \\
{a_{m}}_{3} {b_{p}}_{3} - {a_{m}}_{9} {b_{p}}_{1} & {a_{m}}_{6} {b_{p}}_{3} - {a_{m}}_{9} {b_{p}}_{2} + {b_{p}}_{1} & {b_{p}}_{2} & 0 \\
{b_{p}}_{1} & {b_{p}}_{2} & {b_{p}}_{3} & 0 \\
{b_{m}}_{1} {b_{p}}_{2} - {b_{m}}_{2} {b_{p}}_{1} & {b_{m}}_{1} {b_{p}}_{3} - {b_{m}}_{3} {b_{p}}_{1} & 0 & {b_{p}}_{1}
\end{matrix}\right],
\left[\begin{matrix}
{k_x}_{11} \\ {k_x}_{12} \\ {k_x}_{13} \\ {k_u}_{11}
\end{matrix}\right],
\left[\begin{matrix}
{a_{m}}_{3} - {a_{p}}_{3} \\ {a_{m}}_{6} - {a_{p}}_{6} \\ {a_{m}}_{9} - {a_{p}}_{9} \\ {b_{m}}_{1}
\end{matrix}\right]
\right)
```

