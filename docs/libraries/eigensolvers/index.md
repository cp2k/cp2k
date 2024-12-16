# Eigensolvers

```{toctree}
---
titlesonly:
maxdepth: 2
---
cusolvermp
dlaf
elpa
```

CP2K integrates multiple libraries for the solution of eigenvalue problems.
ScaLAPACK is a mandatory dependency when using MPI, but GPU-accelerated libraries are optionally available:

- [cuSOLVERMp]
- [DLA-Future]
- [ELPA]

[cuSOLVERMp]: https://docs.nvidia.com/cuda/cusolvermp/
[DLA-Future]: https://github.com/eth-cscs/DLA-Future
[ELPA]: https://elpa.mpcdf.mpg.de/
