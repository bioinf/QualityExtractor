#/bin/sh
read -p "Press [Enter] to sort $1 into $1.sorted.bam"
java -jar SortSam.jar I=$1 O=$1.sorted.sam SO=coordinate
./samtools view -bS $1.sorted.sam > $1.sorted.bam
echo "Done."
