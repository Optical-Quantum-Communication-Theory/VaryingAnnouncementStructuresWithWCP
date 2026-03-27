# Performance of BB84 without decoy states under varying announcement structures

This is a public version of the code used in *Performance of BB84 without decoy states under varying announcement structures* \[[arXiv](https://arxiv.org/abs/2603.22448)\]. This was built for [v2.1.0](https://github.com/Optical-Quantum-Communication-Theory/openQKDsecurity/releases/tag/v2.1.0) of the Open QKD Security package.

In the following table we list each main file needed to be run to reproduce the data for each figure:

| Figure  | Main Files Used |
|---------|-----------------|
| Fig. 3  | `BB84_Misalign_Asym.m`, `NPAB_Misalign_Asym.m`, `SARG04_Misalign_Asym.m` |
| Fig. 4  | One can use the asym loss file, fix the loss, and manually change the photon number cutoff.
| Fig. 5  | `BB84_Loss_Asym.m`, `NPAB_Loss_Asym.m`, `SARG04_Loss_Asym.m`|
| Fig. 6  | `BB84_Loss_Asym.m`, `NPAB_Loss_Asym.m`, `SARG04_Loss_Asym.m`|
| Fig. 7  | `BB84_Depol_Asym.m`, `NPAB_Depol_Asym.m`, `SARG04_Depol_Asym.m`|
| Fig. 8  | `BB84_Loss_1e12.m`, `NoPABFiniteLoss1e12.m`, `optimization_SARG_Adaptive_Finite_1e12.m`|
| Fig. 9  | `BB84_Loss_1e9.m`, `NoPABFiniteLoss1e9.m`, `optimization_SARG_Adaptive_Finite_1e9.m`, `optimization_SARG_Adaptive_Finite_1e9_highloss.m`|
| Fig. 10 | `BB84_Loss_1e6.m`, `NoPABFiniteLoss1e6.m`, `optimization_SARG_Adaptive_Finite_1e6.m`|
| Fig. 11 | `BB84_Finite_Depol_1e12.m`, `NoPAB_Depol_1e12.m`, `SARG04_Depol_Finite_1e12.m`|
| Fig. 12 | `BB84_Finite_Depol_1e9.m`, `NoPAB_Depol_1e9.m`, `SARG04_Depol_Finite_1e9.m`|
| Fig. 13 | `BB84_Finite_Misalign_1e12.m`, `NoPAB_Misalign_1e12.m`, `SARG04_Misalign_Finite_1e12.m`|
| Fig. 14 | `BB84_Finite_Misalign_1e9.m`, `NoPAB_Misalign_1e9.m`, `SARG04_Misalign_Finite_1e9.m`|

## Installation instructions
> [!CAUTION]
> This repository is for archival and transparency purposes; we do not guarantee compatibility with other versions of the Open QKD Security package beyond the ones listed above.

This code was designed and tested with Matlab 2024b, though any other recent edition should work. In addition to Matlab, this code requires the following additional resources:
* [Mosek](https://www.mosek.com/) (Optional but it is faster and typically gives better results).
* The Mathworks Statistics and Machine Learning Toolbox.


### As zip
1. Download the linked version of the OpenQKD Security package from above and follow all [installation instructions](https://github.com/Optical-Quantum-Communication-Theory/openQKDsecurity/tree/016911bfe68cd9efb65864b09bd0af0cf3f57c28).
2. Also follow the additional Mosek install instructions if you want an exact match.
3. Download the latest release on the side bar and unzip in your preferred directory and add this folder to the Matlab path.


### with git
1. Clone this repository and its exact submodules navigate to your desired directory and run,
```
git clone --recurse-submodules https://github.com/Optical-Quantum-Communication-Theory/VaryingAnnouncementStructuresWithWCP
```
2. Follow all further [installation instructions](https://github.com/Optical-Quantum-Communication-Theory/openQKDsecurity/tree/016911bfe68cd9efb65864b09bd0af0cf3f57c28).
3. Also follow the additional Mosek install instructions if you want an exact match.
