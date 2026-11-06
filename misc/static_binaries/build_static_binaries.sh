#!/bin/sh
# Builds statically linked bam2cov and hmm_flagger binaries that run on any Linux
# machine of the same architecture (no htslib or other libraries needed at runtime).
#
# Meant to run inside an Alpine container (musl libc makes fully static linking safe):
#   docker run --rm -v "$PWD":/src -w /src alpine:3.20 sh misc/static_binaries/build_static_binaries.sh <out_dir>
#
# The binaries are written to <out_dir> (default: dist) as bam2cov and hmm_flagger.
set -eu

OUT_DIR=${1:-dist}
HTSLIB_VERSION=${HTSLIB_VERSION:-1.13}
PROGRAMS="bam2cov hmm_flagger"

apk add --no-cache build-base git wget perl bash autoconf automake \
    zlib-dev zlib-static bzip2-dev bzip2-static xz-dev xz-static

SRC_DIR=$(pwd)
mkdir -p /home/apps

# sonLib (the programs Makefile expects it in /home/apps/sonLib)
if [ ! -f /home/apps/sonLib/lib/sonLib.a ]; then
    git clone --depth 1 https://github.com/benedictpaten/sonLib /home/apps/sonLib
    make -C /home/apps/sonLib
fi

# htslib as a static library in /usr/local, without libcurl/plugins so it needs no extra shared libraries
if [ ! -f /usr/local/lib/libhts.a ]; then
    cd /tmp
    wget -q https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2
    tar -xjf htslib-${HTSLIB_VERSION}.tar.bz2
    cd htslib-${HTSLIB_VERSION}
    ./configure --disable-libcurl --disable-gcs --disable-s3 --disable-plugins --without-libdeflate
    make -j"$(nproc)" lib-static
    make install
    cd "${SRC_DIR}"
fi

# Build in a copy of programs/ so the checkout is left untouched; remove any objects or binaries
# from earlier builds in the checkout so everything is compiled here against musl
BUILD_DIR=$(mktemp -d)
cp -r programs "${BUILD_DIR}/"
cd "${BUILD_DIR}/programs"
find . -name '*.o' -delete
rm -rf bin
mkdir -p bin
# musl gives new threads only 128 KB of stack by default; with -z stack-size it uses 8 MB
# (the glibc default) instead, so the multi-threaded code behaves as it does with glibc.
STATIC_LINK='-w $(STATIC_LIBS) $(INC) -static -Wl,-z,stack-size=8388608 -lhts -lz -lbz2 -llzma -lm -lpthread'
for prog in ${PROGRAMS}; do
    make -j"$(nproc)" "bin/${prog}" CCFLAGS="${STATIC_LINK}"
done

cd "${SRC_DIR}"
mkdir -p "${OUT_DIR}"
for prog in ${PROGRAMS}; do
    cp "${BUILD_DIR}/programs/bin/${prog}" "${OUT_DIR}/${prog}"
    strip "${OUT_DIR}/${prog}"
    "${OUT_DIR}/${prog}" -h > /dev/null
done
rm -rf "${BUILD_DIR}"
ls -l "${OUT_DIR}"
