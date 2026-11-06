#!/bin/bash
# Builds bam2cov and hmm_flagger for macOS (run it on the architecture you want, e.g. arm64).
#
# macOS does not support fully static binaries, so htslib, liblzma and sonLib are linked in
# statically and the binaries only use libraries that are part of every macOS install
# (libSystem, libz, libbz2). They run on macOS ${MACOSX_DEPLOYMENT_TARGET} or newer.
#
# Needs the Xcode command line tools (xcode-select --install) and curl. Run from the repo root:
#   bash misc/static_binaries/build_macos_binaries.sh <out_dir>
#
# The binaries are written to <out_dir> (default: dist) as bam2cov and hmm_flagger.
set -euo pipefail

OUT_DIR=${1:-dist}
HTSLIB_VERSION=${HTSLIB_VERSION:-1.13}
XZ_VERSION=${XZ_VERSION:-5.6.4}
export MACOSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET:-11.0}
PROGRAMS="bam2cov hmm_flagger"
# Dependencies are built here (reused by later runs)
DEPS_DIR=${DEPS_DIR:-${PWD}/.macos_deps}
NCPU=$(sysctl -n hw.ncpu)

SRC_DIR=${PWD}
PREFIX=${DEPS_DIR}/prefix
mkdir -p "${DEPS_DIR}" "${PREFIX}"

# liblzma as a static library only (it is not part of macOS)
if [ ! -f "${PREFIX}/lib/liblzma.a" ]; then
    cd "${DEPS_DIR}"
    curl -fsSL -o xz-${XZ_VERSION}.tar.gz \
        https://github.com/tukaani-project/xz/releases/download/v${XZ_VERSION}/xz-${XZ_VERSION}.tar.gz
    tar -xzf xz-${XZ_VERSION}.tar.gz
    cd xz-${XZ_VERSION}
    ./configure --prefix="${PREFIX}" --disable-shared --enable-static \
        --disable-xz --disable-xzdec --disable-lzmadec --disable-lzmainfo \
        --disable-lzma-links --disable-scripts --disable-doc --disable-nls
    make -j"${NCPU}"
    make install
fi

# htslib as a static library, without libcurl/plugins so it needs no extra libraries.
# It is used from its build directory, so no libhts.dylib is around to be picked up by the linker.
HTS_DIR=${DEPS_DIR}/htslib-${HTSLIB_VERSION}
if [ ! -f "${HTS_DIR}/libhts.a" ]; then
    cd "${DEPS_DIR}"
    curl -fsSL -o htslib-${HTSLIB_VERSION}.tar.bz2 \
        https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2
    tar -xjf htslib-${HTSLIB_VERSION}.tar.bz2
    cd "${HTS_DIR}"
    ./configure --disable-libcurl --disable-gcs --disable-s3 --disable-plugins --without-libdeflate \
        CPPFLAGS="-I${PREFIX}/include" LDFLAGS="-L${PREFIX}/lib"
    make -j"${NCPU}" lib-static
fi

# sonLib
SONLIB_DIR=${DEPS_DIR}/sonLib
if [ ! -f "${SONLIB_DIR}/lib/sonLib.a" ]; then
    rm -rf "${SONLIB_DIR}"
    git clone --depth 1 https://github.com/benedictpaten/sonLib "${SONLIB_DIR}"
    make -C "${SONLIB_DIR}"
fi

# Build in a copy of programs/ so the checkout is left untouched; remove any objects or binaries
# from earlier builds in the checkout so everything is compiled here
BUILD_DIR=$(mktemp -d)
cd "${SRC_DIR}"
cp -r programs "${BUILD_DIR}/"
cd "${BUILD_DIR}/programs"
find . -name '*.o' -delete
rm -rf bin
mkdir -p bin
# Point the Makefile at the dependencies built above instead of /home/apps/sonLib and /usr/local
sed -i '' \
    -e "s|/home/apps/sonLib|${SONLIB_DIR}|g" \
    -e "s|-I /usr/local/include/htslib|-I ${HTS_DIR} -I ${HTS_DIR}/htslib|" \
    Makefile
# Static libraries are given by path so the linker cannot pick a dylib with the same name
LINK='-w $(STATIC_LIBS) $(INC) '"${HTS_DIR}/libhts.a ${PREFIX}/lib/liblzma.a"' -lz -lbz2 -lm -lpthread'
for prog in ${PROGRAMS}; do
    make -j"${NCPU}" CC=clang "bin/${prog}" CCFLAGS="${LINK}"
done

cd "${SRC_DIR}"
mkdir -p "${OUT_DIR}"
for prog in ${PROGRAMS}; do
    cp "${BUILD_DIR}/programs/bin/${prog}" "${OUT_DIR}/${prog}"
    strip "${OUT_DIR}/${prog}"
    # Apple Silicon only runs signed code; re-sign (ad hoc) after strip
    codesign --force --sign - "${OUT_DIR}/${prog}"
    "${OUT_DIR}/${prog}" -h > /dev/null
done
rm -rf "${BUILD_DIR}"
ls -l "${OUT_DIR}"
