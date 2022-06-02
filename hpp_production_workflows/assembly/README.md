# HPRC Assembly Production Workflow

The HPRC assemblies are produced using Hifiasm which is a fast haplotype-resolved de-novo assembler for PacBio HiFi reads. The trio binning mode of Hifiasm v0.11 is used for this phase of assembly production.

Here is the link to the github of Hifiasm: [github.com/chhylp123/hifiasm](github.com/chhylp123/hifiasm) 

Since the parental short reads are available for each sample, Hifiasm can produce a pair of haplotype-resolved assemblies with trio binning. To perform such an assembly, we counted the k-mers of paternal and maternal reads separately and indexed them using Yak then passed them to Hifiasm for haplotype phasing.

Here is the link to the github of Yak used for k-mer counting and indexing: [github.com/lh3/yak](github.com/lh3/yak)
