FROM debian:jessie
RUN apt-get update && apt-get install -y \
    sudo \
    packaging-dev \
    pkg-config \
    libglib2.0-dev \
    giflib-dbg \
    build-essential \
    automake  && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
RUN git clone --recursive https://github.com/leomrtns/super_sptree.git && \
    mkdir build && \
    cd build && \
    ../super_sptree/configure && \
    make && make install
