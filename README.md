# Lattice Boltzmann reference map technique (LBRMT)

[![arXiv](https://img.shields.io/badge/arXiv-2402.12696-b31b1b?logo=arXiv&logoColor=arXiv&link=https%3A%2F%2Farxiv.org%2Fabs%2F2402.12696)](https://arxiv.org/abs/2402.12696)
![Made with C++](https://img.shields.io/badge/Made_with-C%2B%2B-logoColor?logo=C%2B%2B&logoColor=C%2B%2B)
![GitHub search hit counter](https://img.shields.io/github/search/yue-sun/lbrmt/Hit?logo=github&logoColor=github)


Yue Sun and Chris H. Rycroft, *A fully-integrated lattice Boltzmann method for fluid–structure interaction*, preprint (2024) [[arXiv preprint](https://arxiv.org/abs/2402.12696)]

## What is LBRMT?

The LBRMT is designed to simulate fluid–structure interaction (FSI).

As an Eulerian method, it only requires one fixed computational grid for both solids and fluids. The LBRMT consists of two integral parts: the first involves an Eulerian representation of finite-strain solids, and the second involves a boundary condition for moving deformable interfaces across different densities.
We use the reference map technique (RMT) to model solid large deformation. Additionally, we have developed a new lattice Boltzmann (LB) boundary condition called "smooth flux correction" to model the interfaces between squishy solids and fluids.

The LBRMT is suitable for parallelization and capturing the multi-body contact of hundreds of solids. It can simulate soft solids rotating, settling, floating, and mixing. We believe it could be a promising tool for simulating the spatiotemporal patterns of active matter and studying both the collective and individual dynamics of biological systems.


## Why LBRMT?

- **_One grid to run them all_**: no extra meshes or remeshing needed to track solids
- **_From Boltzmann to lattice Boltzmann_**: easy to parallelize for fluid simulation
- **_The more, the merrier_**: efficient to model hundreds of soft solids


## _"What is a video worth?"_

If a picture is worth a thousand words, check out our [simulation videos](https://github.com/yue-sun/lbrmt/videos/README.md)!


## How to use LBRMT?

Code and tutorial will be released soon.


## More on the RMT

### Software

| Code  | Paper |
| ------------- | ------------- |
| [![IncRMT GitHub](https://img.shields.io/badge/chr1shr_-incrmt-logo?logo=github&link=https%3A%2F%2Fgithub.com%2Fchr1shr%2Fincrmt)](https://github.com/chr1shr/incrmt)  | *Reference map technique for incompressible fluid-structure interaction*  |
| [![Static Badge](https://img.shields.io/badge/ylunalin_-RMT3D-logo?logo=github&link=https%3A%2F%2Fgithub.com%2Fylunalin%2FRMT3D)](https://github.com/ylunalin/RMT3D)  | *Eulerian simulation of complex suspensions and biolocomotion in three dimensions*  |


### References
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

5. Chris H. Rycroft, C.-H. Wu, Y. Yu, and K. Kamrin, *Reference map technique for incompressible fluid-structure interaction*, Journal of Fluid Mechanics **898**, A9 (2020).
   [doi:10.1017/jfm.2020.353](https://doi.org/10.1017/jfm.2020.353)

6. Yuexia Luna Lin, Nicolas J. Derr, and Chris H. Rycroft, *Eulerian simulation of complex suspensions and biolocomotion in three dimensions*, Proc. Natl. Acad. Sci. **119**, e2105338118 (2022). [doi:10.1073/pnas.2105338118](https://doi.org/10.1073/pnas.2105338118)

7. Xiaolin Wang, Ken Kamrin, and Chris H. Rycroft, *An Eulerian method for mixed soft and rigid body interactions in incompressible fluids*, Physics of Fluids **34**, 033604 (2022). [doi:10.1063/5.0082233](https://doi.org/10.1063/5.0082233)
