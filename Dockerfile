FROM continuumio/miniconda3

RUN apt-get -y update \ 
    && apt-get -y install --no-install-recommends \
    procps \
    && rm -rf /var/lib/apt/lists/*

ADD environment.yml /tmp/environment.yml
RUN conda env update -f /tmp/environment.yml
RUN echo "source activate env" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH 