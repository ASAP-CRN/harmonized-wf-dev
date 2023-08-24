FROM ubuntu:jammy AS builder

RUN apt-get -qq update \
	&& apt-get -qq install \
		git \
	&& rm -rf /var/lib/apt/lists/*

ARG CONDA_ENVIRONMENT_REPO_HASH
WORKDIR /opt
RUN git clone https://github.com/shahrozeabbas/Multiome-SeuratV4.git
RUN cd /opt/Multiome-SeuratV4 \
	&& git checkout "${CONDA_ENVIRONMENT_REPO_HASH}"

ARG SCRIPTS_REPO_HASH
RUN git clone https://github.com/shahrozeabbas/Harmony-RNA-Workflow.git
RUN cd /opt/Harmony-RNA-Workflow \
	&& git checkout "${SCRIPTS_REPO_HASH}"


FROM continuumio/miniconda3:23.5.2-0

MAINTAINER Heather Ward <heather@dnastack.com>

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"
ARG CONDA_ENVIRONMENT "${CONDA_ENVIRONMENT}"

# Add conda environment definition
COPY --from=builder /opt/Multiome-SeuratV4/envs/${CONDA_ENVIRONMENT}.yml /opt/envs/
RUN conda env create --file "/opt/envs/${CONDA_ENVIRONMENT}.yml" --name "${CONDA_ENVIRONMENT}"
ENV PATH /opt/conda/envs/${CONDA_ENVIRONMENT}/bin:${PATH}

# Install argparse
RUN Rscript -e 'install.packages("argparse", repos="http://cran.us.r-project.org")'

# Add scripts
COPY --from=builder /opt/Harmony-RNA-Workflow/scripts /opt/scripts

# Add override scripts
COPY scripts/main /opt/scripts/main
COPY scripts/utility /opt/scripts/utility

RUN echo "source activate ${CONDA_ENVIRONMENT}" >> ~/.bashrc
