# Lattice Boltzmann reference map technique (LBRMT)

[![Paper](https://img.shields.io/badge/Paper-10.1016/j.jcp.2025.113774-blue?logoColor=ecf0f1&labelColor=34495e)](https://doi.org/10.1016/j.jcp.2025.113774)
[![arXiv](https://img.shields.io/badge/arXiv-2402.12696-b31b1b?logo=arXiv&logoColor=arXiv&link=https%3A%2F%2Farxiv.org%2Fabs%2F2402.12696)](https://arxiv.org/abs/2402.12696)
![Made with C++](https://img.shields.io/badge/Made_with-C%2B%2B-logoColor?logo=C%2B%2B&logoColor=C%2B%2B)
[![GitHub license](https://img.shields.io/github/license/yue-sun/lbrmt?labelColor=34495e)](https://github.com/yue-sun/lbrmt/blob/main/LICENSE)
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Fyue-sun%2Flbrmt&count_bg=%230D74E7&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=Hits&edge_flat=false)](https://hits.seeyoufarm.com)

Yue Sun and Chris H. Rycroft, *A fully-integrated lattice Boltzmann method for fluid–structure interaction*, Journal of Computational Physics **526**, 113774 (2025) [[DOI link](https://doi.org/10.1016/j.jcp.2025.113774)] [[arXiv preprint](https://arxiv.org/abs/2402.12696)]

<p align="center">
  <img width="640" src=https://github.com/yue-sun/lbrmt/assets/30400317/9546c5ce-d1e1-446e-80d1-ee2c9e68f364>
</p>


## 🌊 What is LBRMT?

The LBRMT is designed to simulate fluid–structure interaction (FSI).

As an Eulerian method, it only requires one fixed computational grid for both solids and fluids. The LBRMT consists of two integral parts: the first involves an Eulerian representation of finite-strain solids, and the second involves a boundary condition for moving deformable interfaces across different densities.
We use the reference map technique (RMT) to model solid large deformation. Additionally, we have developed a new lattice Boltzmann (LB) boundary condition called "smooth flux correction" to model the interfaces between squishy solids and fluids.

The LBRMT is suitable for parallelization and capturing the multi-body contact of hundreds of solids. It can simulate soft solids rotating, settling, floating, and mixing. We believe it could be a promising tool for simulating the spatiotemporal patterns of active matter and studying both the collective and individual dynamics of biological systems.


## 🌟 Why LBRMT?

- **_One grid to run them all_**: no extra meshes or remeshing needed to track solids
- **_From Boltzmann to lattice Boltzmann_**: easy to parallelize for fluid simulation
- **_The more, the merrier_**: efficient to model hundreds of soft solids


## 🎬 _"What is a video worth?"_

If a picture is worth a thousand words, [check out our simulation videos](videos/README.md)!


## 💻 How to use LBRMT?

The code is written in C++ and uses the OpenMP library for multithreading. It has been tested on Linux and MacOS (GCC installed via Homebrew). The following compiling instructions assume that you have a working C++ environment with OpenMP on your computer. All commands are in the terminal.

### ⏬ Install and compile

First, navigate to the directory that you want to install `LBRMT` then clone the repo:
```shell
git clone https://github.com/yue-sun/lbrmt.git
```
Second, navigate into `lbrmt` and compile:
```shell
cd lbrmt
make
```
Now your executables should be in the `build/` and ready to run!
> ⚠️ **GCC**: If you're a Mac user, make sure your `GCC` is in `/opt/homebrew/bin/gcc` when you type `which gcc`. Also, adjust the `GCC` version in the `Makefile` according to your `GCC` version.

### 🆒 Examples

After successfully compiling the code, you should have the following directories within `lbrmt/`:
```shell
lbrmt
│   README.md
│   Makefile
└───build # simulation executables (autogenerated)
│   │   sim_fsi    # exec for lid-driven cavity and settling/floating
│   │   sim_rotate # exec for rotating
│   └───sim_mix    # exec for mixing
└───objs  # objects and linker files (autogenerated)
└───src   # code source files
└───sims  # examples (config files, simulation results also saved here)
│   │   fsi_ldc.cfg     # lid-driven cavity with a soft solid
│   │   fsi_rotate.cfg  # 4 anchored rotors rotating
│   │   fsi_settle.cfg  # 1 solid settling or floating
│   └───fsi_mix.cfg     # 100 solids mixing
└───utils # scripts for visualization and data processing
```

The benchmark example (Section 4.1) of a soft solid in lid-driven cavity can be run using two threads with the following command:
```shell
OMP_NUM_THREADS=2 build/sim_fsi sims/fsi_ldc.cfg
```
The config file specifies the simulation parameters, such as domain size, geometric and material properties for the solid, and physical properties for the fluid. The code will create a directory called `fsi_ldc.out` for the simulation output under the `sims/` directory. To process the output data, navigate to the `utils` directory and there are Python scripts and notebooks to postprocess for data analysis and visualization.

To modify the simulation domain size or solid/fluid properties, modify the `*.cfg` files.

To run other examples, follow the similar command structure:
```shell
OMP_NUM_THREADS=(num_threads) build/(your_exec) sims/(your_sim).cfg
```
> ⚠️ If you run into complication or running errors, please either create a [New issue](https://github.com/yue-sun/lbrmt/issues) on GitHub or contact Yue Sun via email.


## 🚧 Known issues and todos

- [ ] Solids at wall boundaries: Currently can only handle `bctype==3` (fully no-slip box); other wall boundary conditions will be added
- [ ] Consolidate executables
- [ ] Make the entire extrapolation routine parallel
- [ ] Simulate microswimmers


## 📝 Citation

If you use this code in your research or anywhere, please cite the paper:
```bibtex
@article{sun2024fullyintegrated,
  title = {A fully-integrated lattice {{Boltzmann}} method for fluid--structure interaction},
  author = {Sun, Yue and Rycroft, Chris H.},
  year = {2025},
  journal = {Journal of Computational Physics},
  volume = {526},
  pages = {113774},
  issn = {0021-9991},
  publisher = {Elsevier},
  doi = {10.1016/j.jcp.2025.113774}
}
```


## 📑 More on the RMT

### ⌨️ Software

| Code  | Paper |
| ------------- | ------------- |
| [![IncRMT GitHub](https://img.shields.io/badge/chr1shr_-incrmt-logo?logo=github&link=https%3A%2F%2Fgithub.com%2Fchr1shr%2Fincrmt)](https://github.com/chr1shr/incrmt)  | *Reference map technique for incompressible fluid-structure interaction*  |
| [![Static Badge](https://img.shields.io/badge/ylunalin_-RMT3D-logo?logo=github&link=https%3A%2F%2Fgithub.com%2Fylunalin%2FRMT3D)](https://github.com/ylunalin/RMT3D)  | *Eulerian simulation of complex suspensions and biolocomotion in three dimensions*  |

### 📚 References
1. Ken Kamrin, *Stochastic and deterministic models for dense granular flow*,
   Ph.D. thesis, Massachusetts Institute of Technology (2008).
   [DSpace](http://hdl.handle.net/1721.1/43736)

2. Ken Kamrin and Jean-Christophe Nave, *An Eulerian approach to the simulation of deformable solids: application to finite-strain elasticity*,
   [arXiv:0901.3799](https://arxiv.org/abs/0901.3799) (2009).

3. Ken Kamrin, Chris H. Rycroft, and Jean-Christophe Nave, *Reference map
   technique for finite-strain elasticity and fluid–solid interaction*, Journal of the Mechanics and Physics of Solids **60**, 1952–1969 (2012).
   [doi:10.1016/j.jmps.2012.06.003](https://doi.org/10.1016/j.jmps.2012.06.003)

4. Boris Valkov, Chris H. Rycroft, and Ken Kamrin, *Eulerian method for
   multiphase interactions of soft solid bodies in fluids*, Journal of Applied Mechanics **82**, 041011 (2015).
   [doi:10.1115/1.4029765](https://doi.org/10.1115/1.4029765)

5. Chris H. Rycroft, Chen-Hung Wu, Yue Yu, and Ken Kamrin, *Reference map technique for incompressible fluid-structure interaction*, Journal of Fluid Mechanics **898**, A9 (2020).
   [doi:10.1017/jfm.2020.353](https://doi.org/10.1017/jfm.2020.353)

6. Yuexia Luna Lin, Nicolas J. Derr, and Chris H. Rycroft, *Eulerian simulation of complex suspensions and biolocomotion in three dimensions*, Proc. Natl. Acad. Sci. **119**, e2105338118 (2022). [doi:10.1073/pnas.2105338118](https://doi.org/10.1073/pnas.2105338118)

7. Xiaolin Wang, Ken Kamrin, and Chris H. Rycroft, *An Eulerian method for mixed soft and rigid body interactions in incompressible fluids*, Physics of Fluids **34**, 033604 (2022). [doi:10.1063/5.0082233](https://doi.org/10.1063/5.0082233)
