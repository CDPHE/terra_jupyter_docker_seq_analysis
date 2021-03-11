FROM us.gcr.io/broad-dsp-gcr-public/terra-jupyter-base:0.0.19
USER root
#this makes it so pip runs as root, not the user
ENV PIP_USER=false

RUN apt-get update && apt-get install -yq --no-install-recommends \
  python-tk \
  tk-dev \
  libssl-dev \
  xz-utils \
  libhdf5-dev \
  openssl \
  make \
  liblzo2-dev \
  zlib1g-dev \
  libz-dev \
  libcurl4-openssl-dev \
  iqtree \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

ENV HTSLIB_CONFIGURE_OPTIONS="--enable-gcs"

# install samtools
# copied from staphb docker file for samtools
RUN apt-get update && apt-get install -y libncurses5-dev \
  libbz2-dev \
  liblzma-dev \
  libcurl4-gnutls-dev \
  zlib1g-dev \
  libssl-dev \
  gcc \
  wget \
  make \
  perl \
  bzip2
  
RUN mkdir samtools &&\
  mkdir data &&\
  cd samtools &&\
  wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 &&\
  tar -xjf samtools-1.10.tar.bz2 &&\
  rm samtools-1.10.tar.bz2 &&\
  cd samtools-1.10 &&\
  ./configure &&\
  make &&\
  make install

ENV LC_ALL=C
WORKDIR /data

#install iqtree
# copied from stapgh docker file for iqtree
ARG IQTREE2_VER="2.0.3"

#install dependencies
RUN apt-get update && apt-get install -y \
 wget

# download, uncompress iqtree2 tarball; make /data
RUN wget https://github.com/iqtree/iqtree2/releases/download/v${IQTREE2_VER}/iqtree-${IQTREE2_VER}-Linux.tar.gz && \
 tar -xzvf iqtree-${IQTREE2_VER}-Linux.tar.gz && \
 rm -v iqtree-${IQTREE2_VER}-Linux.tar.gz && \
 mkdir /data

# set PATH and locale settings for singularity compatibility
ENV PATH="/iqtree-${IQTREE2_VER}-Linux/bin:${PATH}"\
 LC_ALL=C
WORKDIR /data


RUN pip3 -V \
 && pip3 install --upgrade pip \
 && pip3 install numpy==1.15.2 \
 && pip3 install py4j==0.10.7 \
 && python3 -mpip install matplotlib==3.0.0 \
 && pip3 install pandas==0.25.3 \
 && pip3 install pandas-gbq==0.12.0 \
 && pip3 install pandas-profiling==2.4.0 \
 && pip3 install seaborn==0.9.0 \
 && pip3 install python-lzo==1.12 \
 && pip3 install google-cloud-bigquery==1.23.1 \
 && pip3 install google-api-core==1.6.0 \
 && pip3 install google-cloud-bigquery-datatransfer==0.4.1 \
 && pip3 install google-cloud-datastore==1.10.0 \
 && pip3 install google-cloud-resource-manager==0.30.0 \
 && pip3 install google-cloud-storage==1.23.0 \
 && pip3 install scikit-learn==0.20.0 \
 && pip3 install statsmodels==0.9.0 \
 && pip3 install ggplot==0.11.5 \
 && sed -i 's/pandas.lib/pandas/g' /usr/local/lib/python3.7/dist-packages/ggplot/stats/smoothers.py \
 # the next few `sed` lines are workaround for a ggplot bug. See https://github.com/yhat/ggpy/issues/662
 && sed -i 's/pandas.tslib.Timestamp/pandas.Timestamp/g' /usr/local/lib/python3.7/dist-packages/ggplot/stats/smoothers.py \
 && sed -i 's/pd.tslib.Timestamp/pd.Timestamp/g' /usr/local/lib/python3.7/dist-packages/ggplot/stats/smoothers.py \
 && sed -i 's/pd.tslib.Timestamp/pd.Timestamp/g' /usr/local/lib/python3.7/dist-packages/ggplot/utils.py \
 && pip3 install bokeh==1.0.0 \
 && pip3 install pyfasta==0.5.2 \
 && pip3 install markdown==2.4.1 \
 && pip3 install pdoc3==0.7.2 \
 && pip3 install biopython==1.72 \
 && pip3 install bx-python==0.8.2 \
 && pip3 install fastinterval==0.1.1 \
 && pip3 install matplotlib-venn==0.11.5 \
 && pip3 install bleach==1.5.0 \
 && pip3 install cycler==0.10.0 \
 && pip3 install h5py==2.7.1 \
 && pip3 install html5lib==0.9999999 \
 && pip3 install joblib==0.11 \
 && pip3 install keras==2.1.6 \
 && pip3 install patsy==0.4.1 \
 && pip3 install protobuf==3.7.1 \
 && pip3 install pymc3==3.10.0 \
 && pip3 install pyparsing==2.2.0 \
 && pip3 install Cython \
 && pip3 install pysam==0.15.4 --no-binary pysam \
 && pip3 install python-dateutil==2.6.1 \
 && pip3 install pytz==2017.3 \
 && pip3 install pyvcf==0.6.8 \
 && pip3 install pyyaml==5.3.1 \
 && pip3 install scipy==1.2 \
 && pip3 install tensorflow==2.0.0a0 \
 && pip3 install theano==0.9.0 \
 && pip3 install tqdm==4.19.4 \
 && pip3 install werkzeug==0.12.2 \
 && pip3 install certifi==2017.4.17 \
 && pip3 install intel-openmp==2018.0.0 \
 && pip3 install mkl==2018.0.3 \
 && pip3 install readline==6.2 \
 && pip3 install setuptools==42.0.2 \
 && pip3 install wheel \
 && pip3 install ete3

ENV USER jupyter-user
USER $USER
#we want pip to install into the user's dir when the notebook is running
ENV PIP_USER=true

# Note: this entrypoint is provided for running Jupyter independently of Leonardo.
# When Leonardo deploys this image onto a cluster, the entrypoint is overwritten to enable
# additional setup inside the container before execution.  Jupyter execution occurs when the
# init-actions.sh script uses 'docker exec' to call run-jupyter.sh.
ENTRYPOINT ["/usr/local/bin/jupyter", "notebook"]
