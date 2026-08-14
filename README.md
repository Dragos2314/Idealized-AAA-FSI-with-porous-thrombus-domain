# svMultiPhysics: Inverse_darcy_permeability ignored in FSI equations
Minimal reproducer for a bug in SimVascular/svMultiPhysics: the Inverse_darcy_permeability domain property (Navier-Stokes-Brinkman drag) is applied correctly inside a standalone fluid equation, but has no effect when the same domain is part of a coupled FSI equation. No warning or error is emitted.

The practical consequence: a porous domain in an FSI simulation behaves as free-flowing fluid. In our application (abdominal aortic aneurysm hemodynamics) the intraluminal thrombus is modeled as a Darcy-porous medium, and in FSI its velocities match the free lumen instead of dropping by orders of magnitude.

Suspected root cause

fsi.cpp calls the fluid element routines with the permeability argument hardcoded to 0.0 rather than reading the per-domain property:

File	Line	Call
fsi.cpp	211	fluid::fluid_3d_m(...) with permeability 0.0
fsi.cpp	238	fluid::fluid_2d_m(...) with permeability 0.0

Compare fluid.cpp:556, which performs the correct lookup:

cpp
eq.dmn[cDmn].prop.at(PhysicalProperyType::inverse_darcy_permeability)

A maintainer has confirmed on the community forum that the 0.0 was set intentionally, originally in the context of an unimplemented valve / FSI-contact case, but that the reasoning is undocumented. Forum thread: <ADD LINK>

Version and environment
svMultiPhysics commit 042a1552 (git describe: March_2025-81-g042a1552)
Built and run inside an Apptainer container (svmp-trilinos.sif, Trilinos/MueLu backend), built 2026-06-16
University of Calgary ARC HPC cluster, SLURM, MPI
Both cases in this repo use the built-in fsils linear algebra backend, so a Trilinos-enabled build is not required to reproduce
Contents
mesh_2domain/          lumen (id 1) + porous thrombus (id 3)
mesh_3domain/          lumen (id 1) + wall (id 2) + porous thrombus (id 3)
fluid_only/
  solver_fluid.xml     CASE A (control): fluid equation, porous thrombus
fsi/
  solver_fsi.xml       CASE B (bug):     FSI equation, same porous thrombus

Both meshes are single conforming meshes with shared nodes at all domain interfaces, so no projection between meshes is involved. Domain membership is assigned through domain_ids.dat.

The two cases

Case A and Case B use the same porous domain (ModelRegionID 3) with the same material properties. Domain 3 is byte-for-byte identical between the two XML files:

xml
<Domain id="3" >
   <Equation> fluid </Equation>
   <Density> 1.06e-3 </Density>
   <Viscosity model="Constant" >
      <Value> 4e-3 </Value>
   </Viscosity>
   <Inverse_darcy_permeability> 2.747e8 </Inverse_darcy_permeability>
</Domain>

The inlet flux, timestep size, step count and spectral radius are also identical. The only meaningful difference is <Add_equation type="fluid"> versus <Add_equation type="FSI">, plus the solid wall domain that the FSI case requires.

Both cases are stripped to a minimum: constant inlet flux (parabolic profile), zero-pressure outlet, 20 timesteps. There is no RCR, no Robin support, no prestress, no initialization files and no pulsatile waveform. The runs are not intended to be physically converged - the effect is visible in the velocity field within a few timesteps.

Permeability value K^-1 = 2.747e8 mm^-2 is taken from Adolph et al. 1997 (intraluminal thrombus). Units throughout are mm-g-s: density 1.06e-3 g/mm^3, dynamic viscosity 4e-3 g/(mm s).

Running

Paths inside each XML are relative to the directory containing it, so run each case from inside its own directory:

bash
cd fluid_only
mpirun -n <NPROCS> svmultiphysics solver_fluid.xml

cd ../fsi
mpirun -n <NPROCS> svmultiphysics solver_fsi.xml

With a container:

bash
apptainer exec svmp-trilinos.sif mpirun -n <NPROCS> svmultiphysics solver_fluid.xml

Each case writes VTU output every 5 steps into its own results directory.

Expected versus observed
Case	Equation	Domain 3 permeability	Expected velocity in domain 3	Observed
A	fluid	2.747e8	orders of magnitude below the lumen	as expected - drag applied
B	FSI	2.747e8	orders of magnitude below the lumen	comparable to the lumen - drag not applied

Measured maximum velocity magnitude in domain 3 at step 20:

Case A (fluid): <FILL IN> mm/s
Case B (FSI): <FILL IN> mm/s
Checking the result

In ParaView: open the result VTU, apply a Threshold filter on Domain_ID to isolate the thrombus region, then colour by Velocity magnitude. Compare the range against the lumen region.

Note that svMultiPhysics writes the Domain_ID array only into the first output VTU file of a run; later files in the same run do not carry it. Use the first output file when thresholding, or copy the array across.

A scripted check, printing the maximum velocity magnitude per domain:

python
import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy

reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName("result_001.vtu")   # first output file carries Domain_ID
reader.Update()
grid = reader.GetOutput()

vel = vtk_to_numpy(grid.GetPointData().GetArray("Velocity"))
speed = np.linalg.norm(vel, axis=1)

dom_cell = vtk_to_numpy(grid.GetCellData().GetArray("Domain_ID"))
for d in np.unique(dom_cell):
    ids = vtk.vtkIdList()
    sel = np.where(dom_cell == d)[0]
    pts = set()
    for c in sel:
        grid.GetCellPoints(int(c), ids)
        for k in range(ids.GetNumberOfIds()):
            pts.add(ids.GetId(k))
    pts = np.fromiter(pts, dtype=int)
    print("Domain_ID %d: max |v| = %.4e mm/s" % (d, speed[pts].max()))

Run the same script on the Case A and Case B outputs. In Case A the porous domain shows a large drop relative to the lumen; in Case B it does not.

Background

This came out of a study comparing two hypotheses of abdominal aortic aneurysm development, using idealized geometries in which the intraluminal thrombus is represented as a Darcy-porous medium and the arterial wall as a hyperelastic solid. Modeling thrombus permeability correctly matters for intra-thrombus seepage and for the wall shear stress and oxygen transport quantities derived from it, so an FSI simulation that silently drops the Brinkman drag gives thrombus flow fields that are not physical.
