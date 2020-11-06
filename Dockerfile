FROM continuumio/miniconda3:latest

COPY install_in_docker.sh /root/wg-blimp-setup.sh

RUN /root/wg-blimp-setup.sh
