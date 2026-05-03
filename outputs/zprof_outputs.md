# ZprofOutput: Expected Output Columns

Each output file corresponds to one **phase** × one **vz direction** (`.pvz` or `.nvz`), with naming:
```
<basename>.<NNNNN>.<phase>.<pvz|nvz>.zprof
```
Phases: `phase1`, `phase2`, ..., `whole`.
vz direction: `pvz` ($v_z > 0$), `nvz` ($v_z < 0$).

All `LoadSummedZprofOutputData` columns store the horizontal sum $\sum_{j,i} A (\text{quantity})$,
where $A = (\text{phase mask}) \times (\text{vz mask}) \times \Delta y \Delta x$ is the phase+vz-filtered face area.

---

## Fixed columns (always present)

| Column | Description |
|--------|-------------|
| `k` | z-index |
| `x3v` | Cell-centered z coordinate |
| `area` | $\sum_{j,i} A$ |

---

## From `LoadOutputData` (standard variables: `cons,prim,uov`)

Standard conserved, primitive, and user-output variables from Athena++.

---

## From `LoadZprofOutputData`

Currently **empty** (the example `Pturbz` is commented out).

---

## From `LoadSummedZprofOutputData`

### Hydrodynamics

| Variable | Type | Condition | Description |
|----------|------|-----------|-------------|
| `mZflux` | scalar | `NSCALARS_FEEDBACK > 0` | Metal flux: $\sum A (\rho v_z) Z$ |
| `Pturbz` | scalar | always | Turbulent z-pressure: $\sum A \rho v_z^2$ |
| `gz` | scalar | `SELF_GRAVITY_ENABLED` | Gravitational acceleration: $\sum A g_z$, where $g_z = -\partial\Phi/\partial z$ |
| `rhogz` | scalar | `SELF_GRAVITY_ENABLED` | $\sum A \rho g_z$ |

### MHD

| Variable | Type | Condition | Description |
|----------|------|-----------|-------------|
| `Pmag1,2,3` | vector | `MAGNETIC_FIELDS_ENABLED` | Magnetic pressure per component: $\sum A B_i^2/2$ |
| `vA1,2,3` | vector | `MAGNETIC_FIELDS_ENABLED` | Alfvén speed per component: $\sum A B_i/\sqrt{\rho}$ |

### Energy (`NON_BAROTROPIC_EOS`)

| Variable | Type | Condition | Description |
|----------|------|-----------|-------------|
| `cs` | scalar | `NON_BAROTROPIC_EOS` | Sound speed: $\sum A \sqrt{P/\rho}$ |
| `Ekin1,2,3` | vector | always | Kinetic energy per direction: $\sum A \rho v_i^2/2$ |
| `Ekin_flux1,2,3` | vector | always | Kinetic energy flux per direction: $\sum A (\rho v_z) v_i^2/2$ |
| `Eth_flux` | scalar | `NON_BAROTROPIC_EOS` | Enthalpy flux: $\sum A \frac{\gamma}{\gamma-1} P v_z$ |
| `work` | scalar | `NON_BAROTROPIC_EOS` | Work against pressure gradient: $\sum A (-\mathbf{v} \cdot \nabla P)$ |
| `cool_rate` | scalar | cooling enabled | Cooling rate: $\sum A n_H^2 \Lambda(T)$ |
| `heat_rate` | scalar | cooling enabled | Heating rate: $\sum A n_H \Gamma(T)$ |
| `Egrav` | scalar | `SELF_GRAVITY_ENABLED` | Gravitational energy: $\sum A \rho\Phi$ |
| `Egflux1,2,3` | vector | `SELF_GRAVITY_ENABLED` | Gravitational energy flux per direction: $\sum A (\rho v_i) g_i$ |
| `Sz_Bpress` | scalar | `MAGNETIC_FIELDS_ENABLED` | Poynting flux (pressure term): $\sum A B^2 v_z$ |
| `Sz_Btens` | scalar | `MAGNETIC_FIELDS_ENABLED` | Poynting flux (tension term): $\sum A [-B_z (\mathbf{v}\cdot\mathbf{B})]$ |

### Cosmic Rays (per group `ng = 0..NCRG-1`)

$P_c = e_c/3$, $\hat{B}$ denotes the magnetic field direction, $v_{A,z} \equiv$ `v_adv[ng,2]`.

| Variable | Type | Condition | Description |
|----------|------|-----------|-------------|
| `vAz_mag` | scalar | `CR_ENABLED` + `stream_flag` | $\sum A \lvert v_{A,z}\rvert$ |
| `ng-Fc_st1,2,3` | vector | `CR_ENABLED` + `stream_flag` | Streaming CR flux: $\sum A 4P_c v_{s,i}$ |
| `ng-Fs_B` | scalar | `CR_ENABLED` + `stream_flag` | Streaming CR flux along $\hat{B}$: $\sum A 4P_c v_{s,\hat{B}}$ |
| `ng-Fc_adv1,2,3` | vector | `CR_ENABLED` | Advective CR flux: $\sum A 4P_c v_i$ |
| `ng-Fc_diff1,2,3` | vector | `CR_ENABLED` | Diffusive CR flux magnitude: $\sum A 4P_c \lvert v_{d,i}\rvert$ |
| `ng-Fc_B` | scalar | `CR_ENABLED` | Total CR flux along $\hat{B}$: $\sum A F_{c,\hat{B}}$ |
| `ng-Fa_B` | scalar | `CR_ENABLED` | Advective CR flux along $\hat{B}$: $\sum A 4P_c v_{\hat{B}}$ |
| `ng-Fd_B` | scalar | `CR_ENABLED` | Diffusive CR flux along $\hat{B}$: $\sum A \frac{4}{3}e_c\lvert v_{d,\hat{B}}\rvert$ |
| `ng-Veff1,2,3` | vector | `CR_ENABLED` | Effective CR velocity: $\sum A v_\text{lim} F_{c,i}/(4P_c)$ |
| `ng-kappac` | scalar | `CR_ENABLED` | CR diffusivity: $\sum A \kappa_c = v_\text{lim}/\sigma_\text{diff}$ |
| `ng-Ccr2` | scalar | `CR_ENABLED` | CR effective sound speed² per $\rho$: $\sum A \frac{4}{3}\frac{P_c}{\rho}\frac{v_z + v_{A,z}/2}{v_z + v_{A,z}}$ |
| `ng-rhoCcr2` | scalar | `CR_ENABLED` | $\sum A \frac{4}{3}P_c\frac{v_z + v_{A,z}/2}{v_z + v_{A,z}}$ (for mass-weighted mean) |
| `ng-Ceff2` | scalar | `CR_ENABLED` | Effective sound speed²: $\sum A (C_{cr}^2 + c_s^2)$ |
| `ng-Ceff2-vz2` | scalar | `CR_ENABLED` | Wind acceleration denominator: $\sum A (C_{eff}^2 - v_z^2)$ |
| `ng-GradPc_B` | scalar | `CR_ENABLED` | CR pressure gradient along $\hat{B}$: $\sum A \lvert\nabla_{\hat{B}} P_c\rvert = \sum A \frac{4}{3}e_c\lvert v_{d,\hat{B}}\rvert/\kappa_c$ |
| `rho_ion` | scalar | `CR_ENABLED==2` | Ion density: $\sum A \rho_\text{ion}$ |
| `xion` | scalar | `CR_ENABLED==2` | Ion fraction: $\sum A x_\text{ion}$ |
| `ng-cooling_cr` | scalar | `CR_ENABLED==2` + `losses_flag > 0` | CR collisional cooling: $\sum A (-e_c \lambda_c n_H)$ |
| `ng-heating_cr` | scalar | `CR_ENABLED` + `stream_flag` | CR streaming heating: $\sum A [-\mathbf{v}_s \cdot \boldsymbol{\sigma}(\mathbf{F}_c - \tfrac{4}{3}e_c\mathbf{v})]$ |
| `ng-work_cr` | scalar | `CR_ENABLED` | CR adiabatic work: $\sum A [-v_\parallel\sigma_\parallel(F_{c,\parallel} - \tfrac{4}{3}e_c v_\parallel) + \mathbf{v}_\perp\cdot\nabla P_c]$ |
| `ng-work_cr_direct` | scalar | `CR_ENABLED` | Direct CR work via central-difference gradient: $\sum A (-\mathbf{v} \cdot \nabla P_c)$ |
