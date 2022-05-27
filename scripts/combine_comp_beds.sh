while getopts a:b:o:t:m:h flag
do
    case "${flag}" in
	    h) printf "\nThis command is for producing the final bed output of the FLAGGER pipeline\nFor more information check https://github.com/human-pangenomics/hpp_production_workflows/tree/asset/coverage\n\tRun command:\n\tbash combine_comp_beds.sh \ \n\t -b <BEDS_TAR_GZ> \ \n\t -m <COLORS_TEXT> \ \n\t -t <TRACK_NAME> \ \n\t -o <OUTPUT_BED>\n"; exit 1 ;;
        b) bedsTarGz=${OPTARG};;
        o) outputBed=${OPTARG};;
	t) trackName=${OPTARG};;
	m) mapColorsText=${OPTARG};;
    esac
done

mkdir beds
tar --strip-components 1 -xvzf ${bedsTarGz} --directory beds
         
ERROR_BED=beds/$(ls beds | grep "error")
DUPLICATED_BED=beds/$(ls beds | grep "duplicated")
HAPLOID_BED=beds/$(ls beds | grep "haploid")
COLLAPSED_BED=beds/$(ls beds | grep "collapsed")
UNKNOWN_BED=beds/$(ls beds | grep "unknown")

cat ${ERROR_BED} | awk '{print $0"\tErr"}' > beds/all.bed
cat ${DUPLICATED_BED} | awk '{print $0"\tDup"}' >> beds/all.bed
cat ${HAPLOID_BED} | awk '{print $0"\tHap"}' >> beds/all.bed
cat ${COLLAPSED_BED} | awk '{print $0"\tCol"}' >> beds/all.bed
cat ${UNKNOWN_BED} | awk '{print $0"\tUnk"}' >> beds/all.bed
cat beds/all.bed | sort -k1,1 -k2,2n > beds/all.sorted.bed


echo "track name=\"${trackName}\" visibility=1 itemRgb=\"On\"" > ${outputBed}
awk 'FNR==NR{c[$1]=$2;next}{print $0"\t0\t+\t"$2"\t"$3"\t"c[$4]}' \
    ${mapColorsText} \
    beds/all.sorted.bed \
    >> ${outputBed}


rm -rf beds

