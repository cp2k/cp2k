FROM ubuntu:rolling

# author: Ole Schuett

RUN apt-get update && apt-get install -y --no-install-recommends \
    python                \
    ca-certificates       \
    wget                  \
    unzip                 \
    less                  \
    make                  \
    rsync                 \
    libfindbin-libs-perl  \
&& rm -rf /var/lib/apt/lists/*

WORKDIR /opt/
RUN wget -O cp2k-master.zip https://github.com/cp2k/cp2k/archive/master.zip && \
    unzip cp2k-master.zip  && \
    rm cp2k-master.zip

WORKDIR /opt/cp2k-master/cp2k
#TODO change test_formatting.sh
RUN sed -i s/svn/#svn/g ./tools/formatting/test_formatting.sh
RUN ./tools/formatting/test_formatting.sh

WORKDIR /opt/cp2k_test_formatting
COPY cmd.sh .
CMD ["/opt/cp2k_test_formatting/cmd.sh"]

#EOF