FROM ubuntu:26.04

# author: Ole Schuett

WORKDIR /opt/cp2k-precommit
COPY install_requirements.sh requirements.txt ./
RUN ./install_requirements.sh
ENV PATH="/opt/venv/bin:/opt/cp2k-precommit:$PATH"

ARG REVISION
ENV REVISION=${REVISION}

COPY . ./
CMD ["./precommit_server.py", "--port=8080"]

#EOF
