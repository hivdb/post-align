FROM ubuntu:18.04
ENV LANG C.UTF-8
RUN apt-get -qy update && \
    apt-get -qqy install python3.8-dev python3-pip
ADD requirements.txt /tmp/
RUN python3.8 /usr/bin/pip3 install --no-cache-dir -r /tmp/requirements.txt && rm -rf /tmp/requirements.txt
WORKDIR /app
