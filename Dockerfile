FROM ubuntu:16.04

MAINTAINER Georgy Shapchits <gogi.soft.gm@gmail.com>

RUN mkdir /src

ADD b_system.py src/b_system.py
ADD main.py src/main.py
ADD mpi.py src/mpi.py
ADD requirements.txt src/requirements.txt

RUN apt update -y && \
    apt install -y \
        apt-utils \
        vim \
        iputils-ping \
        build-essential \
        mpich python \
        python-pip \
        python3 \
        python3-pip && \
    pip3 install -r src/requirements.txt

CMD ["sleep", "2d"]
