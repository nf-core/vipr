From:nfcore/base
Bootstrap:docker

%labels
    MAINTAINER Andreas Wilm <wilma@gis.a-star.edu.sg>
    DESCRIPTION Singularity image containing all requirements for the nf-core/vipr pipeline
    VERSION 1.0dev

%environment
    PATH=/opt/conda/envs/nf-core-vipr-1.0dev/bin:$PATH
    export PATH

%files
    environment.yml /

%post
    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a
    pip install git+git://github.com/andreas-wilm/vipr-tools.git@08a360a && pip install git+git://github.com/CSB5/decont.git@bf03c35c
