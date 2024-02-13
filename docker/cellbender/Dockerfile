ARG CELLBENDER_VERSION
FROM us-central1-docker.pkg.dev/dnastack-asap-parkinsons/workflow-images/util:1.1.0 as scripts
FROM us.gcr.io/broad-dsde-methods/cellbender:${CELLBENDER_VERSION}
ENV CELLBENDER_VERSION "${CELLBENDER_VERSION}"

MAINTAINER Karen Fang <karen@dnastack.com>

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG GCLOUD_CLI_VERSION
ENV GCLOUD_CLI_VERSION "${GCLOUD_CLI_VERSION}"
ARG PYTHON_VERSION
ENV PYTHON_VERSION "${PYTHON_VERSION}"
ARG MINICONDA_VERSION
ENV MINICONDA_VERSION "${MINICONDA_VERSION}"
ARG NVIDIA_DRIVER_VERSION
ENV NVIDIA_DRIVER_VERSION "${NVIDIA_DRIVER_VERSION}"
ARG CUDA_VERSION
ENV CUDA_VERSION "${CUDA_VERSION}"

RUN apt-get -qq update && \
	apt-get -qq install \
		time 

COPY ./requirements.txt /opt/requirements.txt
RUN pip install -r /opt/requirements.txt

COPY ./scripts /opt/scripts

COPY --from=scripts /opt/scripts /opt/scripts
ENV PATH "${PATH}:/opt/scripts"
