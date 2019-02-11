FROM debian

ENV PATH="/opt/miniconda3/bin:$PATH"
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

COPY environment.yml /pipeline/environment.yml
COPY install_deps.sh /pipeline/install_deps.sh

SHELL ["/bin/bash", "-c"]

RUN /pipeline/install_deps.sh
