<center> <h1>Codes for a bound-preserving discontinuous Galerkin solver for incompressible two-phase flow problem implemented in Firedrake project</h1> </center>

Numerical examples for
> Mohammad. S. Joshaghani and Beatrice Riviere,
> ``Maximum-principle-satisfying discontinuous Galerkin methods for incompressible two-phase immiscible flow" Available in [arXiv](https://arxiv.org/abs/).

In this repository, we have provided python computer codes for the numerical solution of immiscible two-phase flows in porous media,
which is obtained by augmenting interior penalty DG formulation with post-processing flux and slope limiters. 
The proposed method is accurate, robust, mesh-independent, mass-conservative, and maximum-principle satisfying.
This repo entails several examples of pressure-driven flow and quarter-five spot that account for gravity, capillary effects, 
and heterogeneity. More details are discussed in the paper.

## Notes on limiters
At each time step, following the Newton solver convergence, we first apply a new flux limiter and then a slope limiter. 
Implementation of the flux limiter algorithm is provided in the module ***FluxLimiter*** along with an auxiliary 
flux wrapper module named ***Hsign***. As for the slope limiter, we use the native ***VertexBasedLimiter*** module embedded in the 
Firedrake project.


## Project tree
```
./Codes
├── Pressure_driven
│   ├── Barrier
│   │   ├── FD_Limited.py
│   │   ├── FD_mesh_noLimiter.py
│   │   ├── Limiter
│   │   │   ├── flux_limiter.py
│   │   │   └── hsign.py
│   │   └── mesh
│   │       ├── crack_3.geo
│   │       └── crack_3.msh
│   ├── Gravity
│   │   ├── FD_incompressible.py
│   │   └── Limiter
│   │       ├── flux_limiter.py
│   │       └── hsign.py
│   ├── Homogen
│   │   ├── Compare
│   │   │   ├── FD_incompressible_new.py
│   │   │   └── Limiter
│   │   │       ├── flux_limiter.py
│   │   │       └── hsign.py
│   │   └── NIPGvsUpwind
│   │       ├── FD_incompressible.py
│   │       └── Limiter
│   │           ├── flux_limiter.py
│   │           └── hsign.py
│   └── Nonhomogen
│       ├── FD_incompressible.py
│       └── Limiter
│           ├── flux_limiter.py
│           └── hsign.py
├── Quarter_five_spot
│   ├── Gravity
│   │   ├── FD_well_incompressible.py
│   │   └── Limiter
│   │       ├── flux_limiter_well.py
│   │       └── hsign.py
│   ├── Homogen
│   │   ├── FD_well_incompressible+SL.py
│   │   ├── FD_well_incompressible.py
│   │   ├── FD_well_incompressible_nolimiter.py
│   │   └── Limiter
│   │       ├── flux_limiter_well.py
│   │       └── hsign.py
│   └── SPE10
│       ├── FD_well_incompressible.py
│       └── Limiter
│           ├── flux_limiter_well.py
│           └── hsign.py
└── Verification
    ├── 1D_Buckley_Leverett
    │   ├── BL_Implicit.py
    │   └── Limiter
    │       ├── flux_limiter.py
    │       └── hsign.py
    ├── 2D_BUckley_Leverett
    │   ├── Implicit.py
    │   └── Limiter
    │       ├── flux_limiter.py
    │       └── hsign.py
    └── Convergence
        ├── DG+FL+SL.py
        ├── DG+FL+SL_avg.py
        ├── DG+FL.py
        ├── DG+SL.py
        ├── DG.py
        ├── DG_avg.py
        ├── Limiter
        │   ├── flux_limiter_well.py
        │   └── hsign.py
        ├── batch_DG+FL+SL.sh
        ├── batch_DG+FL+SL_avg.sh
        ├── batch_DG+FL.sh
        ├── batch_DG+SL.sh
        ├── batch_DG.sh
        ├── batch_DG_avg.sh
        └── rate.sh

27 directories, 50 files
```


![GIF](./Video/Video1a.gif)
