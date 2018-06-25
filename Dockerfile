FROM mjmg/centos-supervisor-base:latest

MAINTAINER  Jiucheng <11217007@zju.edu.cn>
ENV REFRESHED_AT=2018-03-20 \
    PIPELINE_VERSION=1.0 \
    SOFT_PATH=/home/software \
    PIPE_PATH=/home/CASS

# env set
RUN yum update -y && \
    yum upgrade -y && \
    yum install -y epel-release \
                   ncurses-devel \
                   ncurses \
                   gcc \
                   g++ \
                   zlib \
                   zlib-devel \
                   bzip2 \
                   xz \
                   xz-devel \
                   bzip2-devel \
                   caja-image-converter.x86_64 \
                   tree \
                   ImageMagick


# install perl 5.16
RUN yum install -y perl \
                   perl-App-cpanminus.noarch  && \
    cpanm Statistics::Basic@1.6611 && \
    cpanm Statistics::Distributions@1.02 && \
    cpanm --force Statistics::Sequences::Runs@0.21 && \
    cpanm --force Statistics::Data@0.11 && \
    cpanm List::MoreUtils@0.416 && \
    cpanm PerlIO::gzip@0.18 && \
    cpanm SVG@2.64

# install openjdk8 and R
RUN yum install -y R && \
    echo "r <- getOption('repos'); r['CRAN'] <- 'https://mirrors.ustc.edu.cn/CRAN/'; options(repos = r);" > ~/.Rprofile


# install conda 
ENV PATH=$SOFT_PATH/bin:$PATH \
    LANG=en_US.UTF-8

COPY ./bssh_native_app-0.9.1.0.tar.gz ./cass_app.py /


RUN wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda.sh && \
    wget https://sourceforge.net/projects/bio-bwa/files/latest/download?source=typ_redirect -O bwa.tar.bz2 && \
    wget https://sourceforge.net/projects/samtools/files/samtools/1.7/samtools-1.7.tar.bz2/download -O samtools.tar.bz2 && \
    sh miniconda.sh -b -f -p $SOFT_PATH && \ 
    conda config --add channels 'https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/' && \
    conda config --set show_channel_urls yes && \
    pip install bssh_native_app-0.9.1.0.tar.gz && \
    pip install sh && \
    pip install simplejson && \
    conda clean -p && conda clean -t && \
    mkdir -p bwa samtools && \
    tar --strip-components=1 -jxvf bwa.tar.bz2 -C bwa && \
    tar --strip-components=1 -jxvf samtools.tar.bz2 -C samtools && \
    make -C bwa && \
    make -C samtools && \
    ln -s /bwa/bwa $SOFT_PATH/bin/bwa && \
    ln -s /samtools/samtools $SOFT_PATH/bin/samtools && \
    rm -rf bssh_native_app-0.9.1.0.tar.gz miniconda.sh  bwa.tar.bz2  samtools.tar.bz2 

COPY ./CASS/* $PIPE_PATH/ 
