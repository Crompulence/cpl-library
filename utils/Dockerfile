# Use latest Ubuntu - focal
FROM ubuntu:20.04
LABEL maintainer="Edward Smith <edward.smith@brunel.ac.uk>"

#Install compilers, mpi (with ssh) and wxpython
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    gcc gfortran build-essential git-core \
    mpich openssh-server \
    python3-dev python3-pip python3-tk python3-matplotlib \ 
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install --upgrade pip \
&&  pip3 install mpi4py numpy pytest scipy

#Clone code from github
RUN git clone https://github.com/Crompulence/cpl-library.git /cpl-library
WORKDIR /cpl-library

#Install CPL library
RUN make PLATFORM=gcc

#Add to the path (same code as SOURCEME.sh)
ENV CPL_PATH=/cpl-library
ENV CPL_BIN_PATH="$CPL_PATH/bin"
ENV PATH=${CPL_BIN_PATH}:$PATH
ENV CPL_LIBRARY_PATH="$CPL_PATH/lib"
ENV CPL_INCLUDE_PATH="$CPL_PATH/include"
ENV LD_LIBRARY_PATH="$CPL_LIBRARY_PATH/:${LD_LIBRARY_PATH}"
ENV PYTHONPATH="$CPL_PATH/src/bindings/python:$PYTHONPATH"
ENV PYTHONPATH="$CPL_PATH/utils:$PYTHONPATH"
ENV CPLPY_PATH="$CPL_PATH/src/bindings/python"

ENV DISPLAY :0
EXPOSE 22

CMD ["bash"]
