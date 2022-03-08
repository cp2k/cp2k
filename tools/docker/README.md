# Docker Containers

This directory hosts docker files for testing and running cp2k.
The only prerequisite is to [install Docker](https://docs.docker.com/get-docker/).

## Production Images

 1. Build the container image:

    <!-- markdownlint-disable MD013 -->
    ```shell
    docker build -f Dockerfile.prod_psmp --build-arg GIT_COMMIT_SHA=$(git rev-parse HEAD) -t cp2k_prod_psmp ../../
    ```
    <!-- markdownlint-enable MD013 -->

 2. Go to the directory with your input files:

    ```shell
    cd ../../benchmarks/QS_single_node
    ```

 3. Run the container:

    <!-- markdownlint-disable MD013 -->
    ```shell
    docker run --shm-size=1g -ti -v "$(pwd)":/mnt cp2k_prod_psmp mpiexec -genv OMP_NUM_THREADS=2 -np 16 cp2k dbcsr.inp
    ```
    <!-- markdownlint-enable MD013 -->

## Test Images

The `Dockerfile.test_*` file are used by the [cp2k-ci](https://github.com/cp2k/cp2k-ci)
to check pull requests and populate the [dashboard](https://dashboard.cp2k.org).

Building a test image executes the test:

```shell
docker build --shm-size=1g -f Dockerfile.test_sdbg -t cp2k_test_sdbg ../../
```

To retrieve the cached report of an old image simply run it:

```shell
docker run cp2k_test_sdbg
```
