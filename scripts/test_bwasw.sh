BWASW="/home/fpavetic/bwa-0.7.4/bwa bwasw"

if [ $# -ne 2 ]
then
    echo "Usage: ./test_bwasw.sh <path-to-bwa-index> <path-to-reads-directory>"
    exit 1
fi

INDEX=$1
INDEX_CHECK="$1.sa"
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

STDOUT="${READS_DIR}/stdout.bwasw"
echo $STDOUT

STDERR="${READS_DIR}/stderr.bwasw"
echo $STDERR

echo "Reads dir: ${READS_DIR}"

for READS in `ls -v $READS_DIR/*.fq`
do
    echo "Processing reads file: $READS"
    RESULT="$READS.bwasw.sam"
    COMMAND="/usr/bin/time -p $BWASW -t 24 -f $RESULT $INDEX $READS"
    echo "Command: $COMMAND"
    $COMMAND
done > $STDOUT 2> $STDERR
