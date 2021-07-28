# versioned base image
FROM rocker/tidyverse:4.0.5

# -- metadata --
LABEL maintainer="Kelsey Montgomery <kelsey.montgomery@sagebase.org>"
LABEL base_image="debian:stretch"
LABEL about.summary="Docker image for sageseqr R package"
LABEL about.license="SPDX:Apache-2.0"

# required to suppress verbose S3 messages
ENV _R_S3_METHOD_REGISTRATION_NOTE_OVERWRITES_=false

# install required packages
RUN apt-get update -y\
&& apt-get install -y dpkg-dev zlib1g-dev libssl-dev libgsl-dev libffi-dev\
&& apt-get install -y curl libcurl4-openssl-dev libglpk-dev libxt6\
&& apt-get install -y git\
&& R -e "install.packages('RcppZiggurat')"\
&& R -e "remotes::install_github('Sage-Bionetworks/sageseqr')"\
&& R -e "install.packages('visNetwork')"
