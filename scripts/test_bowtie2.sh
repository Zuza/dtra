BOWTIE2="~/bowtie2-2.1.0/bowtie2"

if [ $# -ne 2 ]
then
    echo "Usage: ./test_bowtie2.sh <path-to-bowtie2-index> <path-to-reads-directory>"
    exit 1
fi

INDEX=$1
INDEX_CHECK="$1.1.bt2"
if [ $INDEX != "" -a -f $INDEX_CHECK ]
then
  echo "$INDEX found :)"
else
  echo "$INDEX not found :("
  exit 1
fi

READS_DIR=$2
if [ $READS_DIR != "" -a -d $READS_DIR ]
then
  echo "${READS_DIR} found :)"
else
  echo "${READS_DIR} not found :("
  exit 1
fi

STDOUT="${READS_DIR}/stdout.bowtie2"
echo $STDOUT

STDERR="${READS_DIR}/stderr.bowtie2"
echo $STDERR


echo "Reads dir: ${READS_DIR}"

for READS in `ls -v $READS_DIR/*.fq`
do
    echo "Processing reads file: $READS"
    RESULT="$READS.bowtie2.sam"
    COMMAND="/usr/bin/time -p $BOWTIE2 -p 24 -x $INDEX -U $READS -S $RESULT"
    echo "Command: $COMMAND"
    bash -c "$COMMAND"
done > $STDOUT 2> $STDERR
