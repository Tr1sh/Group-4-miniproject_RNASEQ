#! /bin/bash
start=$SECONDS

for i in `cat list.txt`
do
	wget http://h3data.cbio.uct.ac.za/assessments/RNASeq/practice/dataset/${i}_R1.fastq.gz
	wget http://h3data.cbio.uct.ac.za/assessments/RNASeq/practice/dataset/${i}_R2.fastq.gz
	
end=$SECONDS
echo "duration: $((end-start)) seconds."
done

