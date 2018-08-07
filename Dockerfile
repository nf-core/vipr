FROM nfcore/base
MAINTAINER Andreas Wilm <wilma@gis.a-star.edu.sg>
LABEL authors="wilma@gis.a-star.edu.sg" \
    description="Docker image containing all requirements for the nf-core/vipr pipeline"

COPY environment.yml /
RUN conda env update -n root -f /environment.yml && conda clean -a
RUN pip install git+git://github.com/andreas-wilm/vipr-tools.git@08a360a && pip install git+git://github.com/CSB5/decont.git@bf03c35c
