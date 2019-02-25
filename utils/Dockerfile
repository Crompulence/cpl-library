# start from base
FROM ubuntu:16.04
MAINTAINER Edward Smith <edward.smith05@imperial.ac.uk>

#Install compilers, mpi (with ssh) and wxpython
RUN apt-get update && apt-get install -y \
    gcc \
    gfortran \
    git-core \
    build-essential \
    mpich \
    openssh-server \
    python-dev \
    python-pip \
    python-tk \
    python-wxgtk3.0 \
 && rm -rf /var/lib/apt/lists/*

RUN pip install \
    matplotlib==2.2 \
    mpi4py \
    numpy \
    pytest \
    scipy

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
