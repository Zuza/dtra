#yersinia
echo "Yersinia"
time ./test_lisa.sh ~/db/bacteria/yersinia/Yersinia_pestis_a1122.GCA_000222975.1.18.dna.nonchromosomal.fa ../yersinia/ ~/db/bacteria/yersinia/tests/
time ./test_bwasw.sh ~/bwa-0.7.4/yersinia ~/db/bacteria/yersinia/tests/
time ./test_bwamem.sh ~/bwa-0.7.4/yersinia ~/db/bacteria/yersinia/tests/

#banana
echo "Banana"
time ./test_lisa.sh ~/db/plants/banana/banana.fasta ../banana20/ ~/db/plants/banana/tests/
time ./test_bwasw.sh ~/bwa-0.7.4/banana ~/db/plants/banana/tests/
time ./test_bwamem.sh ~/bwa-0.7.4/banana ~/db/plants/banana/tests/

#chicken
echo "Chicken"
time ./test_lisa.sh ~/db/animals/chicken/gallus_gallus.fasta ../chicken20/ ~/db/animals/chicken/tests
time ./test_bwasw.sh ~/bwa-0.7.4/chicken ~/db/animals/chicken/tests/
time ./test_bwamem.sh ~/bwa-0.7.4/chicken ~/db/animals/chicken/tests


