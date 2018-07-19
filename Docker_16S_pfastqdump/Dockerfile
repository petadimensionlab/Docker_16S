#利用するUbuntuのイメージ
FROM ubuntu:14.04

# いろいろインストール
RUN apt-get update -y && apt-get install -y wget tar nano

WORKDIR /

ENV VERSION 2.9.1-1
ENV FILE sratoolkit.$VERSION-ubuntu64
ENV TARFILE $FILE.tar.gz

RUN wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/$VERSION/$TARFILE
RUN tar -xvf $TARFILE && \
    rm $TARFILE && \
    mv $FILE /sratool

ADD vsearch_ubuntu /vsearch_ubuntu
