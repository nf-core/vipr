FROM nfcore/base
MAINTAINER Andreas Wilm <wilma@gis.a-star.edu.sg>
LABEL authors="wilma@gis.a-star.edu.sg" \
    description="Docker image containing all requirements for the nf-core/vipr pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nfcore-vipr-1.0-dev/bin:$PATH
RUN pip install git+git://github.com/andreas-wilm/vipr-tools.git@08a360a && pip install git+git://github.com/CSB5/decont.git@bf03c35c
