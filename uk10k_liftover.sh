java -jar $PICARD CreateSequenceDictionary \
      R=~/Researches/trojan_horses_for_pirnas/original_data/genomes/hg38.fa \
      O=~/Researches/trojan_horses_for_pirna/original_data/genomes/hg38.dict

java -jar $PICARD LiftoverVcf \
     I=~/Researches/trojan_horses_for_pirnas/original_data/polymorphism_data/UK10K/UK10K_COHORT.20160215.sites.chr_added.vcf.gz \
     O=~/Researches/trojan_horses_for_pirnas/original_data/polymorphism_data/UK10K/UK10K_COHORT.20160215.sites.chr_added.lifted.vcf.gz\
     CHAIN=~/Researches/trojan_horses_for_pirnas/original_data/chains/hg19ToHg38.over.chain.gz \
     REJECT=~/Researches/trojan_horses_for_pirnas/original_data/polymorphism_data/UK10K/UK10K_COHORT.20160215.sites.chr_added.rejected.vcf.gz\
     R=~/Researches/trojan_horses_for_pirnas/original_data/genomes/hg38.fa \
     CREATE_INDEX=true