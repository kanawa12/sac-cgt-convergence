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
