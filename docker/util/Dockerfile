FROM gcr.io/google.com/cloudsdktool/google-cloud-cli:444.0.0-slim

RUN python3 -m pip install pandas

COPY scripts /opt/scripts

ENV PATH "${PATH}:/opt/scripts"
