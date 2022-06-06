# syntax=docker/dockerfile:1

ARG VERSION=1.0
FROM python:3.8-slim-buster AS OMACASE

LABEL version=$VERSION
LABEL maintainer="Claire Chung clairechung112@gmail.com"

WORKDIR /omacase-docker

COPY requirements.txt .
RUN python3 -m pip install --upgrade setuptools pip wheel
RUN python3 -m pip install -r requirements.txt

COPY . .
RUN python3 -m pip install .

EXPOSE 8050

CMD [ "python3", "-m", "omacase",  "-m", "web" , "-p", "8050" ]