# i-PI

i-PI is a Python interface for ab initio path integral molecular dynamics simulations. i-PI is
composed of a Python server (i-pi itself, that does not need to be compiled but only requires a
relatively recent version of Python and Numpy) that propagates the (path integral) dynamics of the
nuclei, and of an external code that acts as a client and computes the electronic energy and forces.

The i-PI documentation can be found at <https://ipi-code.org>.

Please cite [](#Kapil2016) and [Kapil2018](https://doi.org/10.1016/j.cpc.2018.09.020) if you use
i-PI with CP2K.

Published work using i-PI with CP2K:

- [J. Phys. Chem. Lett. 2020, 11, 9, 3724–3730](https://doi.org/10.1021/acs.jpclett.0c01025)
- [PNAS 01, 22, 2019 116 (4) 1110-1115](https://doi.org/10.1073/pnas.1815117116)
- [Nature Communications 12, 766 (2021)](https://doi.org/10.1038/s41467-021-20914-0)

## Download i-PI from Github

```none
git clone https://github.com/i-pi/i-pi.git
```

To run i-PI, one need to source the environment file from the i-PI directory

```none
source ${PATH_TO_IPI}/env.sh
```

Then one can run i-PI by using

```none
i-pi input.xml > log &
```

There are many input examples in
[i-pi examples folder](https://github.com/i-pi/i-pi/tree/master/examples), which contains different
methods, e.g. NVE, NVT, NPT, PIMD, REMP, etc.

## Run i-PI with INET socket

To use CP2K as the client code using an **internet** domain socket on the
[HOST](#CP2K_INPUT.MOTION.DRIVER.HOST) address “host_address” and on the
[PORT](#CP2K_INPUT.MOTION.DRIVER.PORT) number “port” the following lines must be added to its input
file:

```none
&MOTION
...
   &DRIVER
      HOST host_address
      PORT port
   &END DRIVER
...
&END MOTION
```

In the `input.xml`, one need to use the same host_address and port to CP2K input.

```none
<ffsocket mode='inet' name='driver'>
  <address>host_address</address>
  <port>port</port>
  <latency>0.01</latency>
  <timeout>5000</timeout>
</ffsocket>
```

## Run i-PI with UNIX socket

If instead a [UNIX](#CP2K_INPUT.MOTION.DRIVER.UNIX) domain socket is required then the following
modification is necessary:

```none
&MOTION
...
   &DRIVER
      HOST host_address
      PORT port
      UNIX
   &END DRIVER
   ...
&END MOTION
```

In the **input.xml**, one to specify the mode=unix to enable the communication via UNIX socket.

```none
<ffsocket mode='unix' name='driver'>
  <address>host_address</address>
  <port>port</port>
  <latency>0.01</latency>
  <timeout>5000</timeout>
</ffsocket>
```

## Run I-PI with CP2K on a supercomputer

To run I-PI with CP2K on a supercomputer is also straightforward, one needs to get the hostname
where the i-PI is executed, and then replace this variable in the CP2K input files. Here is an
example script one can use to run I-PI with CP2K on the Daint machine (CSCS).

```bash
HOST=$(hostname)

source ~/i-pi/env.sh

if [ -e simulation.restart ]; then
   sed -i "s/address>[^<]*</address>$HOST</" simulation.restart
   i-pi simulation.restart  >> log.ipi 2>&1 &
else
   sed -i "s/address>[^<]*</address>$HOST</" input.xml
   i-pi input.xml &> log.ipi &
fi
sleep 5

sed -i "s/HOST.*/HOST $HOST/" cp2k.inp

srun cp2k.psmp -i cp2k.inp -o cp2k.out

wait
```
