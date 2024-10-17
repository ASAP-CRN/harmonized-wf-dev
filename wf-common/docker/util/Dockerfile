FROM gcr.io/google.com/cloudsdktool/google-cloud-cli:444.0.0-slim

RUN apt-get -qq update \
	&& apt-get -qq install \
		xxd

COPY scripts /opt/scripts

ENV PATH "${PATH}:/opt/scripts"
