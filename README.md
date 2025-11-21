# End Simulations

End simulations using rat-pac.

## Installation
- Install Root 6.25 and Geant4 11.0
- Install ratpac-two
- Make sure you `source ratpac.sh` inside ratpac-two repository
- `mkdir build && cd build && make install`

## Running
```
export ENDDATA=$(pwd)/data
source end.sh
end <macro.mac> <rat options>
```
Use vis.mac as an example.

## Usage of docker image
Taking `podman` as an example,
```shell
# under the root directory of this repository
# user 0 --> root user in podman; does not use it when you use `docker`
podman run -it --rm -v ./:/ratpac-setup/ENDSim --user 0 -w /ratpac-setup/ENDSim docker.io/ratpac/ratpac-two
```
After running inside the container, compile under the following way,
```shell
cd /ratpac-setup/ENDSim && mkdir build
cd build && cmake -S ../ -B . -DCMAKE_INSTALL_PREFIX=/ratpac-setup/local/
make && make install
```
After installing the package, source the environment variables,
```shell
. /ratpac-setup/local/bin/end.sh
```

To have X from inside the container,
```shell
podman run -it --rm --user 0 \
  -v ./:/ratpac-setup/ENDSim \
  -v /tmp/.X11-unix:/tmp/.X11-unix \
  -e DISPLAY=$DISPLAY \
  docker.io/ratpac/ratpac-two
```
After entering the container,
```shell
end macros/vis.mac --vis
```
The arg `--vis` or `-g` is required. See details in the
[online documentation](https://ratpac.readthedocs.io/en/latest/users_guide/command_interface.html#id2).
