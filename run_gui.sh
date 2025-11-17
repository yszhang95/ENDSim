podman run -it --rm --user 0 \
  -v ./:/ratpac-setup/ENDSim \
  -v /tmp/.X11-unix:/tmp/.X11-unix \
  -e DISPLAY=$DISPLAY \
  docker.io/ratpac/ratpac-two

# inside docker
# end macros/vis.mac --vis
