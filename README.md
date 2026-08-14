# svMultiPhysics: Inverse_darcy_permeability ignored in FSI equations

Minimal reproducer for a bug in [SimVascular/svMultiPhysics](https://github.com/SimVascular/svMultiPhysics):
the `Inverse_darcy_permeability` domain property (Navier-Stokes-Brinkman drag) is
applied correctly inside a standalone `fluid` equation, but has no effect when the
same domain is part of a coupled `FSI` equation. No warning or error is emitted.

The practical consequence: a porous domain in an FSI simulation behaves as free-flowing
fluid. In our application (abdominal aortic aneurysm hemodynamics) the intraluminal
thrombus is modeled as a Darcy-porous medium, and in FSI its velocities match the
free lumen instead of dropping by orders of magnitude.

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
case, but that the reasoning is undocumented. Forum thread: (https://simtk.org/plugins/phpBB/viewtopicPhpbb.php?f=188&t=26305&p=0&start=0&view=&sid=13249f51ea63f1610756c4aeaa95f4f2)

## Version and environment

- svMultiPhysics commit `042a1552` (`git describe`: `March_2025-81-g042a1552`)
- Built and run inside an Apptainer container (`svmp-trilinos.sif`, Trilinos/MueLu
  backend), built 2026-06-16
- University of Calgary ARC HPC cluster, SLURM, MPI
- Both cases in this repo use the built-in `fsils` linear algebra backend, so a
  Trilinos-enabled build is not required to reproduce

## Contents

```
mesh-idealized-fsi/
  2-domain/            lumen (id 1) + porous thrombus (id 3)
  3-domain/            lumen (id 1) + wall (id 2) + porous thrombus (id 3)
fluid_only/
  solver_fluid.xml     CASE A (control): fluid equation, porous thrombus
fsi/
  solver_fsi.xml       CASE B (bug):     FSI equation, same porous thrombus
```

Both meshes are single conforming meshes with shared nodes at all domain
interfaces, so no projection between meshes is involved. Domain membership is
assigned through `domain_ids.dat`.

## The two cases

Case A and Case B use the same porous domain (`ModelRegionID` 3) with the same
material properties. Domain 3 is byte-for-byte identical between the two XML files:

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
The only meaningful difference is `<Add_equation type="fluid">` versus
`<Add_equation type="FSI">`, plus the solid wall domain that the FSI case requires.

Both cases are stripped to a minimum: constant inlet flux (parabolic profile),
zero-pressure outlet, 20 timesteps. There is no RCR, no Robin support, no
prestress, no initialization files and no pulsatile waveform. The runs are not
intended to be physically converged - the effect is visible in the velocity field
within a few timesteps.

Permeability value `K^-1 = 2.747e8 mm^-2` is taken from Adolph et al. 1997
(intraluminal thrombus). Units throughout are mm-g-s: density 1.06e-3flow fields that are not physical.
