# Porting NEMO4.2.1

Differences between NEMO4.2 and DRAKKAR4.2 should be ported in DRAKKAR4.2.1 taking into accounts differences between NEMO4.2 and NEMO4.2.1

Also NEMO4.2.1 has been run without DCM environment on Adastra

## 

## DRAKKAR modifications in 4.2

| routine | difference | ported as is | differences in NEMO4.2.1 |
|--|--|--|--|
| ICE/icerst.F90 | change in name of restart/output (path) | yes | no |
| ICE/icestp.F90 | change in name of restart/output (path) | yes | no |
| OCE/DIA/diaar5.F90 | halosteric ssh | no | yes |
| OCE/DOM/domain.F90 | change in name of restart (path) | no | yes |
| OCE/DOM/dommsk.F90 | partial steps mask | no | yes |
| OCE/DOM/dtatsd.F90 | modify damping init from TS files| no | yes |
| OCE/ICB/icb_oce.F90 | add iceberg calving restarts ?| yes | no |
| OCE/ICB/icbclv.F90 | add iceberg calving restarts ?| yes | no |
| OCE/ICB/icbini.F90 | add iceberg calving restarts ?| no | yes |
| OCE/ICB/icbrst.F90 | add iceberg calving restarts ?| yes | no |
| OCE/ICB/icbstp.F90 | add iceberg calving restarts ?| no | yes |
| OCE/ICB/icbtrj.F90 | add iceberg calving restarts ?| no | yes |
| OCE/IOM/in_out_manager.F90 | change in name of restart/output (path) | yes | no |
| OCE/IOM/iom.F90 | change in name of restart/output (path) | no | yes |
| OCE/IOM/restart.F90 | change in name of restart/output (path) | no | yes |
| OCE/ISF/isf_oce.F90 | allow for multiple ice shelf input files with different frequency | no | yes |
| OCE/ISF/isfpar.F90 | allow for multiple ice shelf input files with different frequency | yes | no |
| OCE/ISF/isfparmlt.F90 | allow for multiple ice shelf input files with different frequency | no | yes |
| OCE/ISF/isfstp.F90 | allow for multiple ice shelf input files with different frequency | no | yes |
| OCE/LBC/lib_mpp.F90 | change in name of restart/output (path) | no | yes |
| OCE/OBS/diaobs.F90 | change in name of restart/output (path) | no | yes |
| OCE/OBS/obs_profiles_def.F90 | change in name of restart/output (path) | no | yes |
| OCE/OBS/obs_readmdt.F90 | modify reference level | yes | no |
| OCE/OBS/obs_surf_def.F90 | change in name of restart/output (path) | no | yes |
| OCE/OBS/obs_write.F90 | change in name of restart/output (path) | yes | no |
| OCE/SBC/sbcblk.F90 | allow climatological forcings | no | yes |
| OCE/SBC/sbcfwb.F90 | output fwprv ? | no | yes |
| OCE/SBC/sbcrnf.F90 | runoffs from multiple files | no | yes |
| OCE/SBC/sbcssr.F90 | add shapiro filter on SSS restoring | no | yes |
| OCE/TRA/trabbl.F90  | use k- criteria instead of depth criteria | no | yes |
| OCE/USR/usrdef_fmask.F90 | config dependant local enhancement of viscosity | no | yes |
| OCE/ZDF/zdfdrg.F90 | boost of top bottom friction from file ? | no | yes |
| OCE/nemogcm.F90 | timing and error/warnings | no | yes |
| OCE/step_oce.F90 | product diagnostics (uT, vS etc.) | yes | no |
| OCE/stpmlf.F90 | product diagnostics (uT, vS etc.) | no | yes |
| OCE/timing.F90 | allow for longer times | no | yes |
