FROM ubuntu:jammy

RUN apt-get update && \
  DEBIAN_FRONTEND='noninteractive' \
  DEBCONF_NONINTERACTIVE_SEEN='true' \
  apt-get install --yes \
    gfortran \
    build-essential \
    git 

RUN useradd \
  --create-home \
  mag_user

USER mag_user

WORKDIR /home/mag_user

RUN git clone 'https://github.com/geodynamics/mag.git'; cd mag/src; make

ENV PATH="/home/mag_user/mag/src:${PATH}"
