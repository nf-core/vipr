From:nfcore/base
Bootstrap:docker

%labels
    MAINTAINER Andreas Wilm <wilma@gis.a-star.edu.sg>
    DESCRIPTION Singularity image containing all requirements for the nf-core/vipr pipeline

%files
    environment.yml /

%post
    /opt/conda/bin/conda env update -n root -f /environment.yml
/opt/conda/bin/conda clean -a
