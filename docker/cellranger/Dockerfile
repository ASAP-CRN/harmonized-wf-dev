FROM ubuntu:jammy

MAINTAINER Heather Ward <heather@dnastack.com>

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"
ARG CONDA_ENVIRONMENT "${CONDA_ENVIRONMENT}"

RUN apt-get -qq update \
	&& apt-get -qq install \
		wget

ARG CELLRANGER_VERSION
# N.B. this must be passed as a build arg (using the utility script: -b CELLRANGER_DOWNLOAD_QUERY_STRING='xxyy')
#   Enter information here: https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
#   The value of CELLRANGER_DOWNLOAD_QUERY_STRING = everything following the ? in the cURL request URL
ARG CELLRANGER_DOWNLOAD_QUERY_STRING
RUN wget -O cellranger-${CELLRANGER_VERSION}.tar.gz \
	"https://cf.10xgenomics.com/releases/cell-exp/cellranger-${CELLRANGER_VERSION}.tar.gz?${CELLRANGER_DOWNLOAD_QUERY_STRING}" \
	&& tar -zxvf cellranger-${CELLRANGER_VERSION}.tar.gz --directory /opt \
	&& rm cellranger-${CELLRANGER_VERSION}.tar.gz

RUN ln -s /opt/cellranger-${CELLRANGER_VERSION}/cellranger /usr/local/bin/
