#利用するUbuntuのイメージ
FROM ubuntu:16.04

WORKDIR /

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN apt-get update -y && \
    apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates curl git tar nano zlib1g-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.5.4-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

ENV TINI_VERSION v0.16.1
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini

ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD [ "/bin/bash" ]

ENV VERSION 2.9.1-1
ENV FILE sratoolkit.$VERSION-ubuntu64
ENV TARFILE $FILE.tar.gz

RUN wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/$VERSION/$TARFILE
RUN tar -xvf $TARFILE && \
    rm $TARFILE && \
    mv $FILE /sratool


ADD vsearch_ubuntu /vsearch_ubuntu
ADD ./ref.tar.bz2 /

RUN conda config --add channels defaults && \
    conda config --add channels conda-forge && \
    conda config --add channels r && \
    conda config --add channels bioconda && \
    conda install -y conda-forge::openjdk fastqc=0.11.7 trimmomatic==0.36=5 vsearch==2.7 mothur

RUN pip install biopython 




