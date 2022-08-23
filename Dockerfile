# syntax=docker/dockerfile:1

FROM python:3.9

## Install system packages
RUN apt-get update && apt-get install -y \
    mafft \
    ncbi-blast+ \
    bedtools

## Install pysam dependencies
RUN apt-get update && apt-get install -y \
    libncurses5-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev

## Compile MINORg wheel
WORKDIR /mnt/chaelab/rachelle/scripts/minorgpy
ADD dist/minorg-0.2.2.0a0.tar.gz .
WORKDIR minorg-0.2.2.0a0

## Install MINORg offline
## Install MINORg dependencies
## numpy needs to be installed before pysam (from pybedtools) but setup.py doesn't control order so we install it separately first
RUN pip install 'numpy==1.19.5' 'click==8.0.4'
RUN python setup.py install

## Install MINORg
COPY dist/minorg-0.2.2.0a0-py3-none-any.whl minorg.whl
RUN pip install --no-cache minorg

## Run minorg
RUN python -c 'import os; print(os.getcwd()); print(os.listdir())'
ENTRYPOINT ["bash"]
# ENTRYPOINT ["minorg"]
# CMD ["python", "./minorg/main.py"]
# CMD ["python", "-c", "'import os; print(os.getcwd()); print(os.listdir())'"]

# ## Install MINORg from web
# RUN python3 pip install --upgrade \
#     --index-url https://test.pypi.org/simple/ \
#     --extra-index-url https://pypi.org/simple/ minorg
