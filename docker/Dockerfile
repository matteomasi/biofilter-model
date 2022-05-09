FROM quay.io/fenicsproject/stable:latest
USER root
RUN sudo pip3 install -U pip &&\
	sudo pip3 install -U numpy &&\
	sudo pip3 install -U scipy &&\
	sudo pip3 install -U pandas &&\
	sudo pip3 install -U pyyaml &&\
	sudo pip3 install --no-binary=h5py h5py
USER root