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
