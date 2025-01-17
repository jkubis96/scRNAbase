#Dockerfile for sc-rnaseq 
FROM ubuntu:20.04

WORKDIR /app
RUN apt-get update
RUN apt-get update && apt-get install -y locales && rm -rf /var/lib/apt/lists/* \
    && localedef -i en_US -c -f UTF-8 -A /usr/share/locale/locale.alias en_US.UTF-8
ENV LANG en_US.utf8

RUN apt-get update


RUN apt -y install python3.9
RUN apt -y install python3-pip
RUN apt-get update
RUN pip3 install pysam==0.16.0.1
RUN pip3 install biopython==1.78
RUN pip3 install umi_tools==1.0.1
RUN pip3 install gdown

RUN apt-get update

RUN DEBIAN_FRONTEND=noninteractive 
RUN apt-get install -y software-properties-common
RUN apt-get update
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu xenial/'
RUN apt-get update
RUN apt-get install -y libtbb-dev
RUN apt-get -y install r-base=3.6.3-2

RUN apt-get -y install openssl
RUN apt-get -y install libcurl4-openssl-dev
RUN apt-get -y install unzip


RUN apt -y install default-jdk
RUN apt-get install wget
RUN apt-get update

RUN apt-get install -y samtools


RUN apt-get update


RUN gdown 1ndAFxTqHUFjhfBEiFuVs-D1SMKBmhfyI \
	&& dpkg -i rna-star_2.7.3a+dfsg-1build2_amd64.deb \
	&& rm rna-star_2.7.3a+dfsg-1build2_amd64.deb
	


RUN gdown 1nQzT2deG9l0Ho_Nj0splNv9kZIIh-gYv \
	&& dpkg -i fastp_0.20.0+dfsg-1build1_amd64.deb \
	&& rm fastp_0.20.0+dfsg-1build1_amd64.deb


RUN apt-get update && apt-get install -y \
	subread


RUN apt-get -y install r-cran-readr
	
RUN Rscript -e "install.packages(c('stringr', 'dplyr', 'doSNOW', 'foreach', 'doParallel'), repos='https://cran.rstudio.com')"

RUN apt-get update



COPY IndexPip ./
RUN chmod +x IndexPip
COPY CountFeatures ./
RUN chmod +x CountFeatures
ENV PATH="/app:$PATH"


COPY scripts/ ./scripts/
RUN chmod +x ./scripts/genome_prep.R
RUN chmod +x ./scripts/gtf_tool.R

COPY requirements_file/ ./requirements_file/
RUN chmod +x ./requirements_file/Adapters.fa


RUN mkdir -p data
RUN mkdir -p genome

WORKDIR /data


CMD ["/bin/bash"]







