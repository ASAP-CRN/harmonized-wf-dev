FROM us-central1-docker.pkg.dev/dnastack-asap-parkinsons/workflow-images/util:1.1.0 as scripts

# Start from nvidia-docker image with drivers pre-installed to use a GPU
FROM nvcr.io/nvidia/cuda:11.7.1-base-ubuntu18.04

MAINTAINER Karen Fang <karen@dnastack.com>

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ENV DOCKER=true \
    CONDA_AUTO_UPDATE_CONDA=false \
    CONDA_DIR="/opt/conda"
ENV PATH "${CONDA_DIR}/bin:${PATH}"

RUN apt-get -qq update && \
	apt-get -qq install --no-install-recommends \
		time \
		curl \
		ca-certificates \
		sudo \
		git

# Install miniconda
ARG MINICONDA_VERSION
ENV MINICONDA_VERSION "${MINICONDA_VERSION}"
RUN curl -so $HOME/miniconda.sh "https://repo.anaconda.com/miniconda/Miniconda3-py37_${MINICONDA_VERSION}-Linux-x86_64.sh" \
	&& chmod +x $HOME/miniconda.sh \
	&& $HOME/miniconda.sh -b -p $CONDA_DIR \
	&& rm $HOME/miniconda.sh

COPY ./requirements.txt /opt/requirements.txt
RUN pip install -r /opt/requirements.txt

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

# Get compiled crcmod for gsutil
RUN conda install -y -c conda-forge crcmod

# Install cellbender and its dependencies
RUN yes | pip install --no-cache-dir -U "git+https://github.com/broadinstitute/CellBender.git@${IMAGE_TAG}" \
	&& conda clean -yaf \
	&& sudo rm -rf ~/.cache/pip
