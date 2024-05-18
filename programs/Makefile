CC=gcc
INC=-I /home/apps/sonLib/C/inc/ -I /usr/local/include/htslib -I submodules/ptAlignment -I submodules/ptBlock -I submodules/cigar_it -I submodules/common -I submodules/track_reader -I submodules/hmm_utils -I submodules/data_types -I submodules/tpool -I submodules/digamma -I submodules/cJSON -I submodules/hmm -I submodules/chunk -I submodules/bias_detector -I submodules/count_data -I submodules/summary_table
STATIC_LIBS=/home/apps/sonLib/lib/sonLib.a
CCFLAGS= -w $(STATIC_LIBS) $(INC) -lm -lhts -lpthread -lz
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

all: ${OBJS} ${BIN_DIR} $(SRC_BINS) $(TEST_BINS) run_tests

bin/%: src/%.c $(OBJS)
	$(CC) -o $@ $< $(OBJS) $(CCFLAGS)

bin/%: tests/%.c $(OBJS)
	$(CC) -o $@ $< $(OBJS) $(CCFLAGS)



${BIN_DIR}:
	@mkdir -p bin


%.o: %.c %.h
	$(CC) -c -o$@ $< $(CCFLAGS)

run_tests:
	printf "#Test_name\tStatus\n" > tests_status.txt
	./bin/test_summary_table; \
		if [ $$? -eq 0 ]; then printf "API:summary_table\tOK\n" >> tests_status.txt; \
				else printf "API:summary_table\tFAILED\n" >> tests_status.txt; fi	
	./bin/test_track_reader; \
		if [ $$? -eq 0 ]; then printf "API:track_reader\tOK\n" >> tests_status.txt; \
				else printf "API:track_reader\tFAILED\n" >> tests_status.txt; fi
	./bin/test_chunks_creator; \
		if [ $$? -eq 0 ]; then printf "API:chunks_creator\tOK\n" >> tests_status.txt; \
				else printf "API:chunks_creator\tFAILED\n" >> tests_status.txt; fi
	./bin/test_chunk_iterator; \
		if [ $$? -eq 0 ]; then printf "API:chunk_iterator\tOK\n" >> tests_status.txt; \
				else printf "API:chunk_iterator\tFAILED\n" >> tests_status.txt; fi
	./bin/test_common; \
		if [ $$? -eq 0 ]; then printf "API:common\tOK\n" >> tests_status.txt; \
				else printf "API:common\tFAILED\n" >> tests_status.txt; fi
	./bin/test_ptBlock; \
		if [ $$? -eq 0 ]; then printf "API:ptBlock\tOK\n" >> tests_status.txt; \
				else printf "API:ptBlock\tFAILED\n" >> tests_status.txt; fi
	./bin/test_bias_detector; \
		if [ $$? -eq 0 ]; then printf "API:bias_detector\tOK\n" >> tests_status.txt; \
                                else printf "API:bias_detector\tFAILED\n" >> tests_status.txt; fi
	# convert bam to cov with annotations
	./bin/bam2cov --bam tests/test_files/bam2cov/ptBlock_bam2cov_test.bam \
		      --annotationJson tests/test_files/bam2cov/ptBlock_bam2cov_test.json \
		      --output tests/test_files/bam2cov/ptBlock_bam2cov_test.output.cov \
		      --runBiasDetection \
		      --minBiasCoverage 1 \
		      --minBiasLength 1 \
		      --baselineAnnotation annot_1;\
		      cmp tests/test_files/bam2cov/ptBlock_bam2cov_test.output.cov \
		          tests/test_files/bam2cov/ptBlock_bam2cov_test.truth.cov; \
		      return_code=$$?; \
		      if [ $$return_code -eq 0 ]; then printf "PROGRAM:bam2cov_with_uncompressed_cov_file\tOK\n" >> tests_status.txt; \
		      		                  else printf "PROGRAM:bam2cov_with_uncompressed_cov_file\tFAILED\n" >> tests_status.txt; fi ;\
		      if [ $$return_code -eq 0 ]; then printf "Testing bam2cov with uncompressed cov file:" && ./bin/print_OK ; \
		                                  else printf "Testing bam2cov with uncompressed cov file:" && ./bin/print_FAIL ; fi
	# convert bam to cov.gz with annotations
	./bin/bam2cov --bam tests/test_files/bam2cov/ptBlock_bam2cov_test.bam \
		      --annotationJson tests/test_files/bam2cov/ptBlock_bam2cov_test.json \
		      --output tests/test_files/bam2cov/ptBlock_bam2cov_test.output.cov.gz \
		      --runBiasDetection \
		      --minBiasCoverage 1 \
		      --minBiasLength 1 \
		      --baselineAnnotation annot_1; \
		      zcmp tests/test_files/bam2cov/ptBlock_bam2cov_test.output.cov.gz \
		           tests/test_files/bam2cov/ptBlock_bam2cov_test.truth.cov.gz; \
		      return_code=$$?; \
		      if [ $$return_code -eq 0 ]; then printf "PROGRAM:bam2cov_with_compressed_cov_file\tOK\n" >> tests_status.txt; \
                                                  else printf "PROGRAM:bam2cov_with_compressed_cov_file\tFAILED\n" >> tests_status.txt; fi ;\
		      if [ $$return_code -eq 0 ]; then printf "Testing bam2cov with compressed cov file:" && ./bin/print_OK ; \
		                                  else printf "Testing bam2cov with compressed cov file:" && ./bin/print_FAIL ; fi
	# convert bam to bed with annotations
	./bin/bam2cov --bam tests/test_files/bam2cov/ptBlock_bam2cov_test.bam \
		      --annotationJson tests/test_files/bam2cov/ptBlock_bam2cov_test.json \
		      --output tests/test_files/bam2cov/ptBlock_bam2cov_test.output.bed \
		      --runBiasDetection \
		      --minBiasCoverage 1 \
		      --minBiasLength 1 \
		      --baselineAnnotation annot_1; \
		      zcmp tests/test_files/bam2cov/ptBlock_bam2cov_test.output.bed \
		           tests/test_files/bam2cov/ptBlock_bam2cov_test.truth.bed; \
	              return_code=$$?; \
		      if [ $$return_code -eq 0 ]; then printf "PROGRAM:bam2cov_with_ucompressed_bed_file\tOK\n" >> tests_status.txt; \
                                                  else printf "PROGRAM:bam2cov_with_ucompressed_bed_file\tFAILED\n" >> tests_status.txt; fi ;\
		      if [ $$return_code -eq 0 ]; then printf "Testing bam2cov with uncompressed bed file:" && ./bin/print_OK ; \
		                        else printf "Testing bam2cov with uncompressed bed file:" && ./bin/print_FAIL ; fi
	# convert cov to bed with both labels
	./bin/coverage_format_converter -i tests/test_files/coverage_format_converter/test_1_with_labels.cov \
                                        -o tests/test_files/coverage_format_converter/test_1_with_labels.output.bed; \
                                        cmp tests/test_files/coverage_format_converter/test_1_with_labels.output.bed \
                                             tests/test_files/coverage_format_converter/truth/test_1_with_labels.truth.bed; \
					return_code=$$?; \
						if [ $$return_code -eq 0 ]; then printf "PROGRAM:bam2cov_with_ucompressed_bed_file\tOK\n" >> tests_status.txt; \
                                                  else printf "PROGRAM:bam2cov_with_ucompressed_bed_file\tFAILED\n" >> tests_status.txt; fi ;\

                                        if [ $$? -eq 0 ]; then printf "Testing coverage_format_converter for cov -> bed (both labels):" && ./bin/print_OK ; \
                                                          else printf "Testing coverage_format_converter for cov -> bed (both labels):" && ./bin/print_FAIL ; fi
	# convert bed to cov with both labels
	./bin/coverage_format_converter -i tests/test_files/coverage_format_converter/test_1_with_labels.bed \
                                        -o tests/test_files/coverage_format_converter/test_1_with_labels.output.cov \
                                        -f tests/test_files/coverage_format_converter/test_1.fa.fai; \
                                        cmp tests/test_files/coverage_format_converter/test_1_with_labels.output.cov \
                                             tests/test_files/coverage_format_converter/truth/test_1_with_labels.truth.cov; \
					     if [ $$? -eq 0 ]; then printf "Testing coverage_format_converter for bed -> cov (both labels):" && ./bin/print_OK ; \
					                       else printf "Testing coverage_format_converter for bed -> cov (both labels):" && ./bin/print_FAIL ; fi
	# convert cov to bed with only truth labels
	./bin/coverage_format_converter -i tests/test_files/coverage_format_converter/test_1_with_only_truth_labels.cov \
                                        -o tests/test_files/coverage_format_converter/test_1_with_only_truth_labels.output.bed \
                                        -f tests/test_files/coverage_format_converter/test_1.fa.fai; \
                                        cmp tests/test_files/coverage_format_converter/test_1_with_only_truth_labels.output.cov \
                                             tests/test_files/coverage_format_converter/truth/test_1_with_only_truth_labels.truth.cov; \
                                        if [ $$? -eq 0 ]; then printf "Testing coverage_format_converter for cov -> bed (only truth labels):" && ./bin/print_OK ; \
                                                          else printf "Testing coverage_format_converter for cov -> bed (only truth labels):" && ./bin/print_FAIL ; fi
	# convert bed to cov with only truth labels
	./bin/coverage_format_converter -i tests/test_files/coverage_format_converter/test_1_with_only_truth_labels.bed \
		                        -o tests/test_files/coverage_format_converter/test_1_with_only_truth_labels.output.cov \
					-f tests/test_files/coverage_format_converter/test_1.fa.fai; \
					cmp tests/test_files/coverage_format_converter/test_1_with_only_truth_labels.output.cov \
					     tests/test_files/coverage_format_converter/truth/test_1_with_only_truth_labels.truth.cov; \
					if [ $$? -eq 0 ]; then printf "Testing coverage_format_converter for bed -> cov (only truth labels):" && ./bin/print_OK ; \
					                  else printf "Testing coverage_format_converter for bed -> cov (only truth labels):" && ./bin/print_FAIL ; fi
	# convert cov to bedgraph
	./bin/coverage_format_converter -i tests/test_files/coverage_format_converter/test_1_with_labels.cov \
		                        -o tests/test_files/coverage_format_converter/test_1_with_labels.w_20.output.bedgraph \
                                        -f tests/test_files/coverage_format_converter/test_1.fa.fai \
	                                -w 20; \
                                        cmp tests/test_files/coverage_format_converter/test_1_with_labels.w_20.output.bedgraph \
                                            tests/test_files/coverage_format_converter/truth/test_1_with_labels.w_20.truth.bedgraph && \
					cmp tests/test_files/coverage_format_converter/test_1_with_labels.w_20.output.high_mapq.bedgraph \
                                            tests/test_files/coverage_format_converter/truth/test_1_with_labels.w_20.truth.high_mapq.bedgraph && \
					cmp tests/test_files/coverage_format_converter/test_1_with_labels.w_20.output.high_clip.bedgraph \
                                            tests/test_files/coverage_format_converter/truth/test_1_with_labels.w_20.truth.high_clip.bedgraph && \
					cmp tests/test_files/coverage_format_converter/test_1_with_labels.w_20.output.region_idx.bedgraph \
                                            tests/test_files/coverage_format_converter/truth/test_1_with_labels.w_20.truth.region_idx.bedgraph;\
                                        if [ $$? -eq 0 ]; then printf "Testing coverage_format_converter for cov -> bedgraph:" && ./bin/print_OK ; \
                                                          else printf "Testing coverage_format_converter for cov -> bedgraph:" && ./bin/print_FAIL ; fi
	# augment with both truth and prediction labels
	./bin/augment_coverage_by_labels -i tests/test_files/augment_coverage_by_labels/test_1.cov \
		                         -t tests/test_files/augment_coverage_by_labels/test_1.truth_labels.bed \
		                         -p tests/test_files/augment_coverage_by_labels/test_1.prediction_labels.bed \
					 -n 4 \
					 -o tests/test_files/augment_coverage_by_labels/test_1.with_labels.output.cov \
					 -f tests/test_files/augment_coverage_by_labels/test_1.fa.fai ;\
					 cmp tests/test_files/augment_coverage_by_labels/test_1.with_labels.output.cov \
					     tests/test_files/augment_coverage_by_labels/truth/test_1.with_labels.truth.cov; \
					 if [ $$? -eq 0 ]; then printf "Testing augment_coverage_by_labels for both truth and prediction labels:" && ./bin/print_OK ; \
                                                           else printf "Testing augment_coverage_by_labels for both truth and prediction labels:" && ./bin/print_FAIL ; fi
	# augment with only truth labels
	./bin/augment_coverage_by_labels -i tests/test_files/augment_coverage_by_labels/test_1.cov \
                                         -t tests/test_files/augment_coverage_by_labels/test_1.truth_labels.bed \
                                         -n 4 \
                                         -o tests/test_files/augment_coverage_by_labels/test_1.with_only_truth_labels.output.cov \
                                         -f tests/test_files/augment_coverage_by_labels/test_1.fa.fai ;\
	                                 cmp tests/test_files/augment_coverage_by_labels/test_1.with_only_truth_labels.output.cov \
                                             tests/test_files/augment_coverage_by_labels/truth/test_1.with_only_truth_labels.truth.cov; \
                                         if [ $$? -eq 0 ]; then printf "Testing augment_coverage_by_labels for only truth labels:" && ./bin/print_OK ; \
                                                           else printf "Testing augment_coverage_by_labels for only truth labels:" && ./bin/print_FAIL ; fi
	# augment with only prediction labels
	./bin/augment_coverage_by_labels -i tests/test_files/augment_coverage_by_labels/test_1.cov \
                                         -p tests/test_files/augment_coverage_by_labels/test_1.prediction_labels.bed \
                                         -n 4 \
                                         -o tests/test_files/augment_coverage_by_labels/test_1.with_only_prediction_labels.output.cov \
                                         -f tests/test_files/augment_coverage_by_labels/test_1.fa.fai ;\
                                         cmp tests/test_files/augment_coverage_by_labels/test_1.with_only_prediction_labels.output.cov \
                                             tests/test_files/augment_coverage_by_labels/truth/test_1.with_only_prediction_labels.truth.cov; \
                                         if [ $$? -eq 0 ]; then printf "Testing augment_coverage_by_labels for only prediction labels:" && ./bin/print_OK ; \
                                                           else printf "Testing augment_coverage_by_labels for only prediction labels:" && ./bin/print_FAIL ; fi
	@printf "\nTest Summary:\n\n"
	@cat tests_status.txt
