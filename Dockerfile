FROM continuumio/miniconda3:latest

COPY wg-blimp-setup.sh /root/wg-blimp-setup.sh

COPY . /tmp/wg-blimp

RUN /root/wg-blimp-setup.sh
