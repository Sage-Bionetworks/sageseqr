# versioned base image
FROM rocker/tidyverse:4.0.5

# -- metadata --
LABEL maintainer="Kelsey Montgomery <kelsey.montgomery@sagebase.org>"
LABEL base_image="debian:stretch"
LABEL about.summary="Docker image for sageseqr R package"
LABEL about.license="SPDX:Apache-2.0"

# required by synapser#293
ENV CRYPTOGRAPHY_DONT_BUILD_RUST=true

# install required packages
RUN apt-get update -y\
&& apt-get install -y dpkg-dev zlib1g-dev libssl-dev libffi-dev\
&& apt-get install -y curl libcurl4-openssl-dev libglpk-dev libxt6\
&& apt-get install -y git\
&& R -e "remotes::install_github('Sage-Bionetworks/sageseqr')"\
&& R -e "install.packages('visNetwork')"
