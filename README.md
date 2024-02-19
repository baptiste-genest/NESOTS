# Non Euclidean Sliced Optimal Transport Sampling
## Baptiste GENEST, Nicolas COURTY, David COEURJOLLY
### Eurographics 24

Here are the instructions to generate blue noise samples on the sphere and hyperbolic spaces in any dimensions or to sample meshes of genus 0 (see spherical sampling) or genus >= 2 (see local hyperbolic sampling)

### Note that we use CLI11 to provide a clear command line API, you can use --help on any executable to see all available options.
Ref : https://github.com/CLIUtils/CLI11

## Building (Linux/macOS)

We use CMake, hence the usual building commands are required

from root folder :
 ```bash
 mkdir build 
 cd build
 cmake ..
 make -j7
 ```

## Visualization

If you want to see the effects of each parameter or simply to see the algorithm in action you can use the option :
``` bash
--viz 
```

to enable visualization with polyscope : https://polyscope.run/

## Blue Noise on the sphere or on the hyperbolic plane

To generate blue noise points on the sphere in dimension 3, a classical use-case would be
```bash
./spherical_bluenoise --dim 3 --sample_size 1000 
```

The process is the same for the hyperbolic case with the exe.

Note that you can sample any distribution by providing samples with the option --target\_measure nu.pts

## Spherical Mesh Sampling
### CEPS Modification 
We use CEPS (https://github.com/MarkGillespie/CEPS) to maps a genus zero mesh to the sphere,
To generate and output the spherical mesh, you must replace the file `CEPS/src/SphericalUniformization/SphericalUniformization.cpp`
by `this_folder/SphericalUniformization.cpp` before compiling and then launch : 
```bash
CEPS_BUILD_FILE/bin/spherical_uniformize input_mesh.obj
```
to generate the mesh.

To then sample a mesh, use 
```bash
./spherical_mesh_sampling --input_mesh data/bunny.obj --sphere_mesh data/spherical-bunny.obj --sample_size 1000
```

## Local Hyperbolic Mesh Sampling


```bash
./local_hyperbolic_mesh_sampling --input_mesh data/macaca.obj --sample_size 1000
```

## Visualizing results

```bash
./vizu_pc --input_mesh data/macaca.obj --input_points /tmp/out.pts
```

## Projective Plane Sampling

The code here only illustrates visually its use but one could export points clouds in the same way as illustrated in the code of `spherical_bluenoise.cpp` for instance.
