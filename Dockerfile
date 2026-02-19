FROM ubuntu:16.04

# install dependencies from pip3

RUN apt update && \
    apt install -y python3 ncbi-blast+ && \
    apt install -y python3-biopython \
                   python3-pip \
                   wget \
                   unzip

ENV CONDA_PLUGINS_AUTO_ACCEPT_TOS=true

# Install dependencies from conda 
RUN cd /usr/local/ && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py38_4.9.2-Linux-x86_64.sh && \
    bash Miniconda3-py38_4.9.2-Linux-x86_64.sh -b -p /usr/local/miniconda && \
    rm Miniconda3-py38_4.9.2-Linux-x86_64.sh && \
    ln -s /usr/local/miniconda/bin/conda /usr/local/bin/ && \
    conda init bash && \
    /bin/bash -c "source /root/.bashrc" && \
    conda install -c bioconda \
    seqtk=1.3 \
    bowtie2=2.4.1 \
    krakenuniq=0.5.7 \
    kallisto=0.46.0 \
    gmap=2020.10.14 \
    snap-aligner=1.0beta.23 \
    samtools=1.7 \
    openssl=1.1.1i \
    bedtools=2.29.2 \
    bwa=0.7.17 \
    mafft=7.475 \
    tabix=0.2.6 && \
    conda clean -afy
    

# install bcftools (build from source because the rest of this container is too old)
RUN apt-get install -y libhts-dev zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev
RUN wget https://github.com/samtools/bcftools/releases/download/1.23/bcftools-1.23.tar.bz2 && \
    tar -xjf bcftools-1.23.tar.bz2 && \
    cd bcftools-1.23 && \
    apt-get install -y libhts-dev zlib1g-dev libbz2-dev liblzma-dev && \
    make && \
    make install && \
    cd .. && \
    rm -rf bcftools-1.23 bcftools-1.23.tar.bz2

# Install Picard 

# Install PrimerClip 
RUN apt-get install -y libgmp-dev zlib1g-dev
RUN wget -qO- https://github.com/commercialhaskell/stack/releases/download/v2.5.1/stack-2.5.1-linux-x86_64.tar.gz | tar xz && \
    mv stack-2.5.1-linux-x86_64/stack /usr/local/bin/stack && \
    rm -rf stack-2.5.1-linux-x86_64
#RUN /usr/local/bin/stack build 
RUN wget https://github.com/michellejlin/covid_swift_pipeline/releases/download/v0.1-alpha/primerclip-deltest.zip && unzip primerclip-deltest.zip && cd /primerclip-deltest/ && /usr/local/bin/stack build&&  /usr/local/bin/stack install && cd .. 

###########
# ANNOVAR #
###########


# http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz




##########
# JAVA 8 #
##########

# Install OpenJDK-8
RUN apt-get update && \
    apt-get install -y openjdk-8-jdk && \
    apt-get install -y ant && \
    apt-get clean;

# Fix certificate issues
RUN apt-get update && \
    apt-get install ca-certificates-java && \
    apt-get clean && \
    update-ca-certificates -f;

# Setup JAVA_HOME -- useful for docker commandline
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
RUN export JAVA_HOME



CMD ["/bin/bash"]
