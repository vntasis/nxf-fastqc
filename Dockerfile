FROM mambaorg/micromamba:0.8.2
LABEL description="Docker image containing all requirements for the nxf-fastqc pipeline"

# Update the system
RUN apt-get update \
	&& apt-get install -y \
	ksh \
	procps \
	libxt-dev \
	git \
	wget \
	&& apt-get clean
# libxt-dev is required to solve the segfault error caused by cairoVersion() in R

COPY environment.yml .
RUN micromamba install -y -n base -f environment.yml && \
    micromamba clean --all --yes
RUN pip install git+https://github.com/ewels/MultiQC.git

# setup faster for linux
RUN wget -P bin https://github.com/angelovangel/faster/releases/download/v0.1.4/x86_64_linux_faster && \
mv bin/x86_64_linux_faster /usr/local/bin/faster && \
chmod +x /usr/local/bin/faster

# Install DSRC FASTQ compressor
RUN  wget --quiet --no-check-certificate https://github.com/refresh-bio/DSRC/releases/download/v2.0.2/dsrc-linux-x64-static.tar.gz &&\
  tar xzvf dsrc-linux-x64-static.tar.gz &&\
  mv bin/dsrc /usr/local/bin && rm -f dsrc-linux-x64-static.tar.gz && chmod +x /usr/local/bin/dsrc


RUN R -e "install.packages('sparkline', repos='http://cran.rstudio.com/')"

# Create a docker user and set home as the working directory
RUN useradd -ms /bin/bash docker
WORKDIR /home/docker
USER docker
ENV BASH_ENV ''

# Run shell
CMD ["bash"]
