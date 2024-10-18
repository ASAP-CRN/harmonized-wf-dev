ARG CELLBENDER_VERSION
FROM us-central1-docker.pkg.dev/dnastack-asap-parkinsons/workflow-images/util:1.1.1 as scripts
FROM us.gcr.io/broad-dsde-methods/cellbender:${CELLBENDER_VERSION}
ENV CELLBENDER_VERSION "${CELLBENDER_VERSION}"

LABEL MAINTAINER="Karen Fang <karen@dnastack.com>"

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
ARG CUDA_VERSION
ENV CUDA_VERSION "${CUDA_VERSION}"

RUN apt-get -qq update && \
	apt-get -qq install \
		time

COPY --from=scripts /opt/scripts /opt/scripts
ENV PATH "${PATH}:/opt/scripts"
