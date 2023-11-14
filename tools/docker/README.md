# Docker Containers for Testing

> **Note** For production containers please visit: <http://github.com/cp2k/cp2k-containers>.

This directory hosts docker files for testing cp2k. They are mostly used by the
[cp2k-ci](https://github.com/cp2k/cp2k-ci) to check pull requests and populate the
[dashboard](https://dashboard.cp2k.org).

To run a test one simply has to build the image:

```shell
docker build --shm-size=1g -f Dockerfile.test_sdbg -t cp2k_test_sdbg ../../
```

To retrieve the cached report of an old image simply run it:

```shell
docker run cp2k_test_sdbg
```
