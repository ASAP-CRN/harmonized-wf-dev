ARG CUDA_VERSION
FROM us-central1-docker.pkg.dev/dnastack-asap-parkinsons/workflow-images/util:1.1.0 as scripts

FROM nvcr.io/nvidia/cuda:${CUDA_VERSION}-base-ubuntu20.04
ENV CUDA_VERSION "${CUDA_VERSION}"

MAINTAINER Karen Fang <karen@dnastack.com>

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"
ARG CUDA_VERSION
ENV CUDA_VERSION "${CUDA_VERSION}"

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -qq update \
	&& apt-get -qq install \
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
		liblzma-dev

ARG PYTHON3_VERSION
ENV PYTHON3_VERSION "${PYTHON3_VERSION}"
RUN curl -O https://www.python.org/ftp/python/${PYTHON3_VERSION}/Python-${PYTHON3_VERSION}.tar.xz && \
	tar -xvf Python-${PYTHON3_VERSION}.tar.xz --directory /opt/ && \
	rm Python-${PYTHON3_VERSION}.tar.xz
RUN cd /opt/Python-${PYTHON3_VERSION} && \
	./configure && \
	make && \
	make altinstall

ENV PATH "${PATH}:/opt/Python-${PYTHON3_VERSION}"

RUN ln -s /opt/Python-${PYTHON3_VERSION}/python /opt/Python-${PYTHON3_VERSION}/python3

COPY ./requirements.txt /opt/requirements.txt
RUN python3 -m pip install -r /opt/requirements.txt

# CUDA-enabled jaxlib is needed
RUN python3 -m pip install jaxlib==0.4.25+cuda12.cudnn89 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

# Add override scripts
COPY scripts/main /opt/scripts/main
COPY scripts/utility /opt/scripts/utility

# Add resources
COPY resources /opt/resources

# gcloud sdk; needed to upload output files
ARG GCLOUD_CLI_VERSION
ENV GCLOUD_CLI_VERSION "${GCLOUD_CLI_VERSION}"
RUN wget "https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-${GCLOUD_CLI_VERSION}-linux-x86_64.tar.gz" \
	&& tar -zxvf "google-cloud-cli-${GCLOUD_CLI_VERSION}-linux-x86_64.tar.gz" --directory /opt \
	&& rm "google-cloud-cli-${GCLOUD_CLI_VERSION}-linux-x86_64.tar.gz"

ENV PATH "${PATH}:/opt/google-cloud-sdk/bin"

COPY --from=scripts /opt/scripts /opt/scripts
ENV PATH "${PATH}:/opt/scripts"
