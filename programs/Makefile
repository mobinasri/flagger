CC=gcc
INC=-I /home/apps/sonLib/C/inc/ -I /usr/local/include/htslib -I submodules/ptAlignment -I submodules/ptBlock -I submodules/cigar_it -I submodules/common -I submodules/block_it -I submodules/hmm -I submodules/data_types -I submodules/tpool -I submodules/digamma -I submodules/cjson
STATIC_LIBS=/home/apps/sonLib/lib/sonLib.a
CCFLAGS=$(STATIC_LIBS) $(INC) -lm -lhts -lpthread -lz
BIN_DIR=bin
SRC_DIR=src
TEST_DIR=tests
MDL_DIR=submodules

MDLS=$(wildcard ${MDL_DIR}/*)
OBJS=$(join $(MDLS), $(MDLS:${MDL_DIR}/%=/%.o))
SRCS=$(wildcard ${SRC_DIR}/*.c ${TEST_DIR}/*.c)
SRC_BINS=$(patsubst ${SRC_DIR}/%.c,${BIN_DIR}/%,$(SRCS))
TEST_BINS=$(patsubst ${TEST_DIR}/%.c,${BIN_DIR}/%,$(SRCS))

.PHONY: all

all: ${OBJS} ${BIN_DIR} $(SRC_BINS) $(TEST_BINS)

bin/%: src/%.c $(OBJS)
	$(CC) -o $@ $< $(OBJS) $(CCFLAGS)

bin/%: tests/%.c $(OBJS)
	$(CC) -o $@ $< $(OBJS) $(CCFLAGS)



${BIN_DIR}:
	@mkdir -p bin


%.o: %.c %.h
	$(CC) -c -o$@ $< $(CCFLAGS)
