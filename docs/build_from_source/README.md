## Building HMM-Flagger from source

This guide is for users who cannot use the `mobinasri/flagger` docker image and want to build the
HMM-Flagger programs themselves on Linux.

The steps are the same on amd64 (x86_64) and arm64 (aarch64) machines. You can check which one you have
with `uname -m`.

### 1. Install system dependencies

You need a C compiler, `make`, autotools and the compression/network libraries that htslib links against.
On Ubuntu/Debian (with root access):

```
sudo apt-get update
sudo apt-get install -y build-essential git wget autoconf automake perl \
    zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev \
    python3 python3-dev python3-pip
```

On other distributions install the equivalent packages (e.g. `gcc make zlib-devel bzip2-devel xz-devel
libcurl-devel openssl-devel` on RHEL/Rocky). Without root access you can get the same libraries from a
conda/mamba environment (`conda install -c conda-forge gcc make zlib bzip2 xz libcurl openssl`).

### 2. Choose an install location

Pick any directory you can write to:

```
export APPS_DIR=${HOME}/flagger_apps      # where sonLib and htslib go
mkdir -p ${APPS_DIR}
```

### 3. Build sonLib and htslib

The HMM-Flagger C programs link against [sonLib](https://github.com/benedictpaten/sonLib) and
[htslib 1.13](https://github.com/samtools/htslib).

```
# sonLib
cd ${APPS_DIR}
git clone https://github.com/benedictpaten/sonLib
cd sonLib
make

# htslib 1.13
cd ${APPS_DIR}
wget https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2
tar -xvjf htslib-1.13.tar.bz2
cd htslib-1.13
autoreconf -i
./configure --prefix=${APPS_DIR}/htslib
make
make install
```

If you have root access, you can skip `--prefix` so that htslib is installed in `/usr/local` and clone
sonLib into `/home/apps/sonLib`; then `programs/Makefile` will find both without any of the environment
variables in the next section.

### 4. Build the HMM-Flagger programs

[`programs/Makefile`](../../programs/Makefile) expects sonLib in `/home/apps/sonLib` and htslib in
`/usr/local`. If you installed them somewhere else (section 3), point the compiler to them with
environment variables and override `STATIC_LIBS` on the `make` command line:

```
# Clone this repository (skip if you already have it)
git clone https://github.com/mobinasri/flagger.git
cd flagger/programs

# Tell gcc where the sonLib/htslib headers and libraries are
export CPATH=${APPS_DIR}/sonLib/C/inc:${APPS_DIR}/htslib/include/htslib:${APPS_DIR}/htslib/include:${CPATH}
export LIBRARY_PATH=${APPS_DIR}/htslib/lib:${LIBRARY_PATH}
export LD_LIBRARY_PATH=${APPS_DIR}/htslib/lib:${LD_LIBRARY_PATH}

make STATIC_LIBS=${APPS_DIR}/sonLib/lib/sonLib.a
```

`make` compiles every program into `programs/bin` and then runs the unit tests; a summary is printed at
the end and saved in `programs/tests_status.txt`. All lines should say `OK`.

### 5. Set up your environment

Add these lines to your `~/.bashrc` so they persist (`LD_LIBRARY_PATH` must include htslib whenever you
run the programs):

```
export FLAGGER_DIR=/path/to/flagger
export PATH=${FLAGGER_DIR}/programs/bin:${PATH}
export LD_LIBRARY_PATH=${APPS_DIR}/htslib/lib:${LD_LIBRARY_PATH}
```

The python scripts in `programs/src` need these packages:

```
pip3 install numpy scipy pandas matplotlib pysam biopython whatshap
# only needed for tune_alpha_hmm_flagger.py
pip3 install smt
```

Check that everything works:

```
hmm_flagger -h
bam2cov -h
```

### Running the python scripts

The [main README](../../README.md) runs the python scripts from `/home/programs/src`, which only exists
inside the docker image. Without docker, use `${FLAGGER_DIR}/programs/src` instead, for example:

```
python3 ${FLAGGER_DIR}/programs/src/filter_hmm_flagger_calls.py \
    --inputBam self_hom_alignment.bam \
    --inputBed final_flagger_prediction.bed \
    --outputBed final_flagger_prediction.conservative.bed \
    --threads 8
```

The self-alignment BAM used here (step 4.1 of the Quick Start) is made with minimap2 and samtools, which
you need to install separately.
