FROM ubuntu:20.04
ENV LANG=C.UTF-8 DEBIAN_FRONTEND=noninteractive
RUN apt-get -qy update && \
    apt-get -qqy install python3.9-dev python3.9-distutils binutils gcc
ADD https://bootstrap.pypa.io/get-pip.py /tmp/get-pip.py
RUN python3.9 /tmp/get-pip.py
COPY requirements.txt /tmp/
RUN /usr/local/bin/pip install --no-cache-dir -r /tmp/requirements.txt && rm -rf /tmp/requirements.txt
COPY /dockerscripts/build.sh /build.sh
RUN chmod +x /build.sh
WORKDIR /app
