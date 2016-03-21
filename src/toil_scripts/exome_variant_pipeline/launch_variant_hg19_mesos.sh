#!/usr/bin/env bash
# John Vivian
#
# Please read the associated README.md before attempting to use.
export PYTHONPATH=$(python -c 'from os.path import abspath as a, dirname as d;import sys;print d(d(d(a(sys.argv[1]))))' $0)
python -m toil_scripts.exome_variant_pipeline.exome_variant_pipeline \
aws:us-west-2:jvivian-releases-exome-1 \
--retryCount 1 \
--config /home/mesosbox/shared/config.txt \
--reference 'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_hg19/hg19.fa' \
--phase 'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_hg19/1000G_phase1.indels.hg19.sites.vcf' \
--mills 'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf' \
--dbsnp 'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_hg19/dbsnp_138.hg19.vcf' \
--cosmic 'https://s3-us-west-2.amazonaws.com/cgl-pipeline-inputs/variant_hg19/cosmic.hg19.vcf' \
--ssec /home/mesosbox/shared/master.key \
--s3_dir cgl-driver-projects/test/variants/ \
--workDir /var/lib/toil \
--batchSystem="mesos" \
--mesosMaster mesos-master:5050 \
#--restart
