# OpenRDP docker container
# Bradley R. Jones
#
# To build:
# docker build -t openrdp .
#
# To run:
# cat alignment.fasta > docker run -i --rm openrdp > output.csv
#
# You may need to use sudo

FROM python:3.8-slim

WORKDIR /app

# install dependencies
RUN apt-get update && \
apt-get install -y build-essential gfortran git wget

# install python dependencies
RUN pip3 install cython numpy scipy matplotlib

# 3seq 
RUN mkdir /git ; \ 
cd /git && \
git clone https://gitlab.com/lamhm/3seq && \
cd 3seq && \
git checkout b854293dec119fd4b56eb17fa1648b478d1233fd && \
make

# geneconv
RUN mkdir /scratch ; \
cd /scratch && \
wget https://www.math.wustl.edu/~sawyer/geneconv/unix.source.tar.gz && \
tar -xzf unix.source.tar.gz && \
cd unix.source && \
gcc -DUNIX -o geneconv -O3 geneconv.c version.c vcalc.c vtcalc.c \
  vsetopts.c vread.c vdump.c vutil.c -lm
 
# copy
RUN mkdir /git/OpenRDP
COPY . /git/OpenRDP/

# install
RUN cd /git/OpenRDP && \
cp /git/3seq/3seq openrdp/bin/3Seq/3seq.Unix && \
cp /scratch/unix.source/geneconv openrdp/bin/GENECONV/geneconv.Unix && \ 
python3 setup.py install && \
ln -s /git/OpenRDP/docker/docker.cfg /app/cfg.ini 

# driver script
# modify docker/docker.cfg to your liking
# add "-m ..." to restrict methods used
RUN echo 'cat - > /app/in.fa && openrdp /app/in.fa --outfile=/app/out.csv --cfg=/app/cfg.ini 1>&2 && cat /app/out.csv' > /app/driver.sh

ENTRYPOINT ["bash", "/app/driver.sh"]

