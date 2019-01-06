FROM alpine:latest

ENV PATH /opt/conda/bin:$PATH

RUN apk --update add --no-cache \
    curl \
    glib \
    libstdc++ \
    nano \
 && curl -L "https://github.com/sgerrand/alpine-pkg-glibc/releases/download/2.27-r0/glibc-2.27-r0.apk" -o /tmp/glibc.apk \
 && apk add --no-cache --allow-untrusted /tmp/glibc.apk \
 && rm -rf /tmp/glibc.apk.sh \
 && curl -L "https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh" -o /tmp/miniconda.sh \
 && sh /tmp/miniconda.sh -b -p /opt/conda \
 && rm /tmp/miniconda.sh \
 && conda upgrade -y pip \
 && conda config --add channels conda-forge \
 && conda config --add channels bioconda \
 && conda install -y conda-forge::openjdk fastqc=0.11.7 trimmomatic==0.36=5 vsearch==2.7 mothur=1.39.5 \
 && conda clean --all \
 && pip install biopython

ADD vsearch_ubuntu /vsearch_ubuntu
ADD ./ref.tar.bz2 /

ENTRYPOINT ["/bin/sh"]