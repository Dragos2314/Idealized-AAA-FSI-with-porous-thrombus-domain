# svMultiPhysics: Inverse_darcy_permeability ignored in FSI equations

Minimal reproducer for a bug in [SimVascular/svMultiPhysics](https://github.com/SimVascular/svMultiPhysics):
the `Inverse_darcy_permeability` domain property (Navier-Stokes-Brinkman drag) is
applied correctly inside a standalone `fluid` equation, but has no effect when the
same domain is part of a coupled `FSI` equation. No warning or error is emitted.

Measured on the two cases in this repository, at the same timestep, with the same
porous domain and the same permeability value:

| Case | Equation | Lumen mean speed | Thrombus mean speed | Ratio |
| --- | --- | --- | --- | --- |
| A | `fluid` | 4.5298e+01 mm/s | 8.4649e-04 mm/s | 53513 |
| B | `FSI` | 3.0638e+01 mm/s | 2.1966e+01 mm/s | 1.39 |

The porous domain damps the flow by four to five orders of magnitude under a fluid
equation, and by essentially nothing under an FSI equation. Thrombus mean velocity
is roughly 26000 times too high in the FSI case.

The practical consequence: a porous domain in an FSI simulation behaves as
free-flowing fluid. In our application (abdominal aortic aneurysm hemodynamics) the
intraluminal thrombus is modeled as a Darcy-porous medium, and in FSI its velocities
match the free lumen instead of dropping.

## Suspected root cause

`fsi.cpp` calls the fluid element routines with the permeability argument hardcoded
to `0.0` rather than reading the per-domain property:

| File | Line | Call |
| --- | --- | --- |
| `fsi.cpp` | 211 | `fluid::fluid_3d_m(...)` with permeability `0.0` |
| `fsi.cpp` | 238 | `fluid::fluid_2d_m(...)` with permeability `0.0` |

Compare `fluid.cpp:556`, which performs the correct lookup:

```cpp
eq.dmn[cDmn].prop.at(PhysicalProperyType::inverse_darcy_permeability)
```

A maintainer has confirmed on the community forum that the `0.0` was set
intentionally, originally in the context of an unimplemented valve / FSI-contact
case, but that the reasoning is undocumented. Forum thread: <ADD LINK>

## Version and environment

- svMultiPhysics commit `042a1552` (`git describe`: `March_2025-81-g042a1552`)
- Built and run inside an Apptainer container (`svmp-trilinos.sif`, Trilinos/MueLu
  backend), built 2026-06-16
- University of Calgary ARC HPC cluster, SLURM, Open MPI, 24 ranks
- Both cases use the built-in `fsils` linear algebra backend, so a Trilinos-enabled
  build is not required to reproduce

## Contents

```
mesh-idealized-fsi/
  2-domain/            lumen (id 1) + porous thrombus (id 3)
  3-domain/            lumen (id 1) + wall (id 2) + porous thrombus (id 3)
fluid_only/
  solver_fluid.xml     CASE A (control): fluid equation, porous thrombus
fsi/
  solver_fsi.xml       CASE B (bug):     FSI equation, same porous thrombus
check_velocity.py      per-domain velocity magnitude summary
```

Both meshes are single conforming meshes with shared nodes at all domain
interfaces, so no projection between meshes is involved. Domain membership is
assigned through `domain_ids.dat`.

## The two cases

Case A and Case B use the same porous domain (`ModelRegionID` 3) with the same
material properties. Domain 3 is identical between the two XML files:

```xml
<Domain id="3" >
   <Equation> fluid </Equation>
   <Density> 1.06e-3 </Density>
   <Viscosity model="Constant" >
      <Value> 4e-3 </Value>
   </Viscosity>
   <Inverse_darcy_permeability> 2.747e8 </Inverse_darcy_permeability>
</Domain>
```

The inlet flux, timestep size, step count and spectral radius are also identical.
The differences are `<Add_equation type="fluid">` versus `<Add_equation type="FSI">`,
the solid wall domain that the FSI case requires, and the mesh-motion equation that
comes with it.

Both cases are stripped to a minimum: constant inlet flux (parabolic profile,
-25286.9 mm^3/s), zero-pressure outlet, 20 timesteps of 1e-3 s. There is no RCR, no
Robin support, no prestress, no initialization files and no pulsatile waveform. The
runs are not intended to be physically converged - the effect is visible in the
velocity field within a few timesteps.

Permeability value `K^-1 = 2.747e8 mm^-2` is taken from Adolph et al. 1997
(intraluminal thrombus). Units throughout are mm-g-s: density 1.06e-3 g/mm^3,
dynamic viscosity 4e-3 g/(mm s).

### Caveat on the comparison

The two cases cannot use an identical mesh, because an FSI equation requires a solid
domain that a rigid-wall fluid case does not have. The 3-domain mesh is the 2-domain
mesh with a wall layer added; the lumen and thrombus geometry and the thrombus
material properties are unchanged between them. The comparison is therefore between
the closest pair of cases the solver permits, not between two runs on one mesh.

## Running

Paths inside each XML are relative to the directory containing it, so run each case
from inside its own directory. With a container:

```bash
apptainer exec --bind <REPO_ROOT>:/case --pwd /case/fluid_only <CONTAINER>.sif \
  mpiexec -n 24 <PATH_TO>/svmultiphysics solver_fluid.xml 2>&1 | tee run_fluid.log

apptainer exec --bind <REPO_ROOT>:/case --pwd /case/fsi <CONTAINER>.sif \
  mpiexec -n 24 <PATH_TO>/svmultiphysics solver_fsi.xml 2>&1 | tee run_fsi.log
```

Without a container:

```bash
cd fluid_only && mpiexec -n 24 svmultiphysics solver_fluid.xml
cd ../fsi     && mpiexec -n 24 svmultiphysics solver_fsi.xml
```

Each case writes VTU output every 5 steps into a `24-procs/` directory. On 24 ranks
the fluid case takes about 33 s/step and the FSI case about 120 s/step.

## Checking the result

svMultiPhysics writes the `Domain_ID` array only into the first output VTU of a run.
With the save stride of 5 used here, that file is `result_005.vtu`, and both cases
are compared at that same step.

Note that the `Domain_ID` values written to the VTU are not the `ModelRegionID`
values set in the XML. The observed mapping is a power-of-two encoding:

| `ModelRegionID` in XML | `Domain_ID` in VTU | Region |
| --- | --- | --- |
| 1 | 2 | lumen |
| 2 | 4 | wall (Case B only) |
| 3 | 8 | thrombus |

Read the unique values out of the file rather than assuming; `check_velocity.py`
does this automatically.

```bash
cd fluid_only/24-procs && python ../../check_velocity.py result_005.vtu
cd ../../fsi/24-procs  && python ../../check_velocity.py result_005.vtu
```

In ParaView the equivalent check is a **Threshold** filter on `Domain_ID` isolating
the value 8, coloured by `Velocity` magnitude.

## Full measured output

Case A, `fluid_only/24-procs/result_005.vtu`:

```
Domain_ID 2: max |v| = 1.6162e+02  mean |v| = 4.5298e+01     (lumen)
Domain_ID 8: max |v| = 8.6793e-03  mean |v| = 8.4649e-04     (thrombus)
```

Case B, `fsi/24-procs/result_005.vtu`:

```
Domain_ID 2: max |v| = 1.9590e+02  mean |v| = 3.0638e+01     (lumen)
Domain_ID 4: max |v| = 3.6589e+01  mean |v| = 9.4639e+00     (wall)
Domain_ID 8: max |v| = 1.4949e+02  mean |v| = 2.1966e+01     (thrombus)
```

Maximum values are reported alongside means because in a conforming mesh the nodes
on the lumen-thrombus interface are shared between both domains and therefore carry
near-lumen velocities in either case. The means are the cleaner discriminator, but
both tell the same story: in Case A the porous domain damps by a factor of 5.4e4,
in Case B by a factor of 1.4.

The Case B run shown here was stopped by a wall-clock limit after 18 of 20
timesteps. This does not affect the measurement, which is taken at step 5.

## Background

This came out of a study comparing two hypotheses of abdominal aortic aneurysm
development, using idealized geometries in which the intraluminal thrombus is
represented as a Darcy-porous medium and the arterial wall as a hyperelastic solid.
Modeling thrombus permeability correctly matters for intra-thrombus seepage and for
the wall shear stress and transport quantities derived from it, so an FSI simulation
that silently drops the Brinkman drag gives thrombus flow fields that are not
physical.
