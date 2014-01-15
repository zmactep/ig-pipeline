for i in {1..35}
do
  for ((j = 10; j <= 1000; j += 10))
  do
  	echo $i " -  " $j
    python ../diff_info.py ../../../../data/test/germline/test-$i/train-$j/prediction.kabat ../../../../data/test/germline/test-$i/test-data/vh-test.kabat > ../../../../data/test/germline/test-$i/train-$j/diff_info.txt
    python ../compare_marking.py ../../../../data/test/germline/test-$i/train-$j/prediction.kabat ../../../../data/test/germline/test-$i/test-data/vh-test.kabat > ../../../../data/test/germline/test-$i/train-$j/compare_marking.txt
  done
done
