FROM continuumio/anaconda3
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN git clone --recursive git@github.com:leomrtns/super_sptree.git
RUN mkdir build
RUN cd build
RUN ../super_sptree-master/configure
RUN make; make install

