#Download base image ubuntu
FROM ubuntu:18.04
run apt-get update
run apt-get -y upgrade
run apt-get -y install --no-install-recommends sudo


###
# Install some basics.
###
run apt-get -y install --no-install-recommends bzip2 \
                                                    g++ \
                                                    gfortran \
                                                    libtool \
                                                    automake \
                                                    autoconf \
                                                    m4 \
                                                    bison \
                                                    flex \
                                                    libcurl4-openssl-dev \
                                                    zlib1g-dev \
                                                    git \
                                                    wget \
                                                    curl \
                                                    libjpeg-dev \
                                                    cmake \
                                                    python-dev \
                                                    cython \
                                                    python-numpy \
                                                    gdb \
                                                    dos2unix \
                                                    antlr \
                                                    libantlr-dev \
                                                    libexpat1-dev \
                                                    libxml2-dev \
                                                    gsl-bin \
                                                    libgsl0-dev \
                                                    udunits-bin \
                                                    libudunits2-0 \
                                                    libudunits2-dev \
                                                    clang \
                                                    zip \
                                                    valgrind \
                                                    python-setuptools \
                                                    make \
                                                    build-essential \
                                                    less \
                                                    unzip \
                                                    patch \
                                                    libsz2 \
                                                    libaec-dev\
                                                    gcc \
                                                    ssh \
                                                    open-coarrays-bin

RUN apt-get -yq install gcc \
                        build-essential \
                        wget \
                        bzip2 \
                        tar \
                        libghc-zlib-dev \
                        m4 \
                        libopenmpi-dev \
                        file \
                        openmpi-bin
RUN apt-get -yq install libcoarrays-dev

# install
RUN apt-get -y install libnetcdf-dev libnetcdff-dev libhdf5-serial-dev  &&\
    apt-get -y install libkernlib1-gfortran netcdf-bin hdf5-tools mpich nano

RUN apt-get -yq install mlocate
RUN updatedb

RUN useradd -ms /bin/bash newuser
WORKDIR /home/newuser
USER newuser







