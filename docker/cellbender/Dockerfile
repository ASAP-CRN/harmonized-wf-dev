FROM us-central1-docker.pkg.dev/dnastack-asap-parkinsons/workflow-images/util:1.1.0 as scripts

# Start from nvidia-docker image with drivers pre-installed to use a GPU
FROM nvcr.io/nvidia/cuda:11.7.1-base-ubuntu18.04

MAINTAINER Karen Fang <karen@dnastack.com>

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -qq update && \
	apt-get -qq install \
		build-essential \
		wget \
		time \
		curl \
		zlib1g-dev \
		libncursesw5-dev \
		libssl-dev \
		libsqlite3-dev \
		tk-dev \
		libgdbm-dev \
		libc6-dev \
		libbz2-dev \
		libffi-dev \
		liblzma-dev \
		ca-certificates \
		python3-dev \
		python3-pip

COPY ./requirements.txt /opt/requirements.txt
RUN python3 -m pip install -r /opt/requirements.txt

# gcloud sdk; needed to upload output files
ARG GCLOUD_CLI_VERSION
ENV GCLOUD_CLI_VERSION "${GCLOUD_CLI_VERSION}"
RUN wget "https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-${GCLOUD_CLI_VERSION}-linux-x86_64.tar.gz" \
	&& tar -zxvf "google-cloud-cli-${GCLOUD_CLI_VERSION}-linux-x86_64.tar.gz" --directory /opt \
	&& rm "google-cloud-cli-${GCLOUD_CLI_VERSION}-linux-x86_64.tar.gz"

ENV PATH "${PATH}:/opt/google-cloud-sdk/bin"

COPY ./scripts /opt/scripts

COPY --from=scripts /opt/scripts /opt/scripts
ENV PATH "${PATH}:/opt/scripts"
