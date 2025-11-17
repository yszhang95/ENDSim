podman run -it --rm -v ./:/ratpac-setup/ENDSim --user 0 docker.io/ratpac/ratpac-two
# # inside docker
# cd /ratpac-setup/ENDSim
# mkdir build
# cd build && cmake -S ../ -B . -DCMAKE_INSTALL_PREFIX=/ratpac-setup/local/
# make && make install
