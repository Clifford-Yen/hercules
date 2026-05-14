# Hercules

## Project History & Maintenance
Hercules was originally developed by researchers at Carnegie Mellon University (CMU). The original source code and theoretical foundations were established by the CMU team. 

Since 2023, this repository has been maintained, modernized, and documented by Clifford (Chu-Han) Yen under the supervision of Prof. Ertugrul Taciroglu at UCLA CEE. The current maintenance efforts focus on code modernization, providing comprehensive documentation, and ensuring compatibility with modern computing environments.

## Authors and Acknowledgments
Please refer to the `AUTHORS.txt` file for the complete list of the original CMU developers. For issues related to the current repository state or documentation, please contact the current maintainer.

## Quick Start
For instructions on how to build, install and execute Hercules, please read
the documentation: https://clifford-yen.github.io/HerculesDoc/index.html 

Some general description of the core components in this directory are given as follows.

| File           | Description |
| -------------- | ----------- |
| `octor`        | Hercules mesh generator/partitioner/load balancer module
| `quake`        | Hercules solver module
| `user.mk`      | User-specified compiler and linker settings. This should be the only file we need to change.
| `systemdef.mk` | generic system-specific compile and linker settings. Tailored by entries in user.mk
| `common.mk`    | system-independent settings
