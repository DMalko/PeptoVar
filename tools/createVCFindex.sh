# create compressed VCF with index
cp ../testdata/test.vcf ../testdata/test.vcf.bak
bgzip -f ../testdata/test.vcf
tabix -p vcf ../testdata/test.vcf.gz
mv ../testdata/test.vcf.bak ../testdata/test.vcf
