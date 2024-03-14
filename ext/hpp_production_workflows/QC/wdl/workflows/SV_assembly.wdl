version 1.0

import "../tasks/parliament.wdl" as runParliament
import "../tasks/sniffles.wdl" as runSniffles
import "../tasks/filterSV.wdl" as runFilterSV
import "../tasks/iris.wdl" as runIris
import "../tasks/jasmine.wdl" as runJasmine

workflow sv_assembly{
    
    meta{
        author: "Avani Khadilkar"
        email: "akhadilk@ucsc.edu"
        description: "WDL implementation of the small variant like errors calling section of the [T2T Polishing Case Study](https://github.com/arangrhie/T2T-Polish/blob/master/doc/T2T_polishing_case_study.md)."
    }
    
    input{
        File IlluminaInputBam
        File IlluminaIndexBam
        File HifiInputBam
        File OntInputBam
        File RefGenome
        File IndexGenome
        String SampleName
        
        Boolean? ParlFilterShortContigs = true
        Boolean? RunParl = false
        String? ParlOtherArgs
        Int? maxDist
        Float? minSeqID
        Int? specReads
    }

    # Run PARLIAMENT on Illumina data if user has chosen to

    if (RunParl == true){
        call runParliament.Parliament as Parl{
            input:
                inputBam = IlluminaInputBam,
                refGenome = RefGenome,
                indexBam = IlluminaIndexBam,
                indexGenome = IndexGenome,
                SampleName = SampleName,
                filterShortContigs = ParlFilterShortContigs,
                otherArgs = ParlOtherArgs
        }
    }

    # Run SNIFFLES on HIFI data

    call runSniffles.Sniffles as HiFiSniffles{
        input:
            inputBam = HifiInputBam,
            SampleName = SampleName,
            outputFileTag = "HiFi"
    }
    
    # Run SNIFFLES on ONT data

    call runSniffles.Sniffles as OntSniffles{
        input:
            inputBam = OntInputBam,
            SampleName = SampleName,
            outputFileTag = "Ont"
    }

    # Run filter.py on HiFi Sniffles Output

    call runFilterSV.Filter as HiFiFilter{
        input:
            inputVcf = HiFiSniffles.vcfOut,
            SampleName = SampleName,
            outputFileTag = "HiFi"
    }

    # Run filter.py on Ont Sniffles Output

    call runFilterSV.Filter as OntFilter{
        input:
            inputVcf = OntSniffles.vcfOut,
            SampleName = SampleName,
            outputFileTag = "Ont"
    }

    call runIris.Iris as Iris{
        input:
            genomeIn = RefGenome,
            readsIn = HifiInputBam,
            vcfIn = HiFiFilter.vcfOut,
            SampleName = SampleName
    }

    call runJasmine.Jasmine as Jasmine{
        input:
            InputVCFs = select_all([Iris.vcfOut, OntFilter.vcfOut, Parl.vcfOut]),
            SampleName = SampleName,
            maxDist = maxDist,
            minSeqID = minSeqID,
            specReads = specReads
    }

    output{
        File? ParlOutput = Parl.vcfOut
        File SnifflesHiFiOutput = HiFiSniffles.vcfOut
        File SnifflesOntOutput = OntSniffles.vcfOut
        File FilterHiFiOutput = HiFiFilter.vcfOut
        File FilterOntOutput = OntFilter.vcfOut
        File IrisHiFiOutput = Iris.vcfOut
        File SV_filelist = Jasmine.SV_filelist
        File SV_like_errors = Jasmine.vcfOut
    }
}
