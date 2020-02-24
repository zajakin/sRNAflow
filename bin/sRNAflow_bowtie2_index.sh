threads=4 
if [ $2 != "" ]
then
	threads=$2 
fi
bowtie2-build --threads $threads $1 `echo $1 | sed 's/.fa$//'`
