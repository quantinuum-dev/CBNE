FROM ubuntu:22.04
LABEL Description="Build environment"

ENV HOME=/root

SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get -y --no-install-recommends install \
    build-essential \
    cmake \
    gdb \
    wget

# Let us add some heavy dependency
RUN cd ${HOME} && \
    wget --no-check-certificate  \
        https://sourceforge.net/projects/boost/files/boost/1.77.0/boost_1_77_0.tar.gz && \
        tar xzf ./boost_1_77_0.tar.gz && \
        cd ./boost_1_77_0 && \
        ./bootstrap.sh && \
        ./b2 install && \
        cd .. && \
        rm -rf ./boost_1_77_0

RUN mkdir /benchmarks

WORKDIR /CBNE

COPY . /CBNE

WORKDIR /CBNE/build/

RUN cmake ..
RUN make

ENTRYPOINT ["/CBNE/build/cbne"]
