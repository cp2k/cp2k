FROM ubuntu:20.04

# author: Ole Schuett

WORKDIR /opt/cp2k-precommit
COPY . /opt/cp2k-precommit/
RUN ./install_requirements.sh

ARG REVISION
ENV REVISION=${REVISION}

CMD ["gunicorn", "--bind=:8080", "--workers=1", "--threads=8", "--timeout=0", "precommit_server:app"]

#EOF
