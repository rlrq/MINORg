# syntax=docker/dockerfile:1

FROM python:3.9

## Install system packages
RUN apt-get update && apt-get install -y \
    mafft \
    ncbi-blast+ \
    bedtools

## Install pysam system dependencies
RUN apt-get update && apt-get install -y \
    libncurses5-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev

## Install MINORg dependencies
## numpy & click need to be installed before pysam (from pybedtools) and typer (respectively) but setup.py doesn't control order so we install it separately first
RUN pip install 'numpy==1.19.5' 'click==8.0.4'

## install the rest here. setup.py gets modified whenever I update the version number, and if we leave the installation of these packages there pybedtools especially will take forever
RUN pip install 'biopython==1.79' 'pybedtools==0.9.0' 'pyfaidx==0.6.4' 'regex' 'typer==0.4.0' 'multiprocess==0.70.12.2'

## Copy files
WORKDIR /minorg_docker
ADD examples_for_docker/examples.zip .
ADD dist/minorg-0.2.3.3a0.tar.gz .
WORKDIR minorg-0.2.3.3a0

## Install MINORg dependencies
RUN python setup.py install

## Install MINORg
COPY dist/minorg-0.2.3.3a0-py3-none-any.whl minorg.whl
RUN pip install --no-cache minorg

## Run minorg
RUN python -c 'import os; print(os.getcwd()); print(os.listdir())'
ENV MINORG_IN_DOCKER Yes
WORKDIR ..
ENTRYPOINT ["bash"]
# ENTRYPOINT ["minorg"]
# CMD ["python", "./minorg/main.py"]
# CMD ["python", "-c", "'import os; print(os.getcwd()); print(os.listdir())'"]

# ## Install MINORg from web
# RUN python3 pip install --upgrade \
#     --index-url https://test.pypi.org/simple/ \
#     --extra-index-url https://pypi.org/simple/ minorg
