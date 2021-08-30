FROM r-base:3.6.3
#FROM rocker/tidyverse:3.6.3

# install requirements
RUN R -e "install.packages('MCMCpack',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('truncnorm',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN apt-get update && apt-get install -y git build-essential libgtkmm-3.0-dev libboost-all-dev libgsl-dev libz-dev libbz2-dev libgsl-dev libcurl4-gnutls-dev cmake
RUN R -e "install.packages('Rfast',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('R.utils',dependencies=TRUE, repos='http://cran.rstudio.com/')"


# copy scripts
COPY ./R/main_docker_multiomic.R /bin/main_docker_multiomic.R
RUN chmod 777 /bin/main_docker_multiomic.R
COPY ./R/function_bulk_repulsive_github.R /bin/function_bulk_repulsive_github.R
RUN chmod 777 /bin/function_bulk_repulsive_github.R

# copy toy data
RUN mkdir /bin/test_data
COPY ./test_data/RNA_dummy.tsv /bin/test_data/RNA_dummy.tsv
COPY ./test_data/proteo_dummy.tsv /bin/test_data/proteo_dummy.tsv
COPY ./test_data/LM22_combined_cell_types.tsv /bin/test_data/LM22_combined_cell_types.tsv

# VOLUME ['/tmp']

WORKDIR '/bin'
