# AMD64 needed explicitly on ARM as some software is only available in AMD64
FROM --platform=amd64 mambaorg/micromamba:2-ubuntu24.04

WORKDIR /home/mambauser
# Copy backend and environment.yml
COPY --chown=mambauser:mambauser environment.yml .
COPY --chown=mambauser:mambauser requirements.txt .

USER mambauser
ENV HOME=/home/mambauser
ENV ENV_NAME=mcquac

RUN echo 'show_banner: false' > ~/.mambarc

RUN micromamba env create -y -f environment.yml \
    && micromamba clean --all --yes


USER root
# First is necessary for base_image to actvate the conda environment second is entrypoint
# which adds the python file to PATH
ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh"]
