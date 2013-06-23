LISA=~/dtra/bin/client
LISA_OPTS="--validate_wgsim=report --long_read_algorithm=lis"

if [ $# -ne 3 ]
then
    echo "Usage: ./test_lisa.sh <path-to-fasta-db> <path-to-lisa-index> <path-to-reads-directory>"
    exit 1
fi

DB=$1
if [ $DB != "" -a -f $DB ]
then
  echo "$DB found :)"
else
  echo "$DB not found :("
  exit 1
fi

INDEX=$2
if [ $INDEX != "" -a -d $INDEX ]
then
  echo "$INDEX found :)"
else
  echo "$INDEX not found :("
  exit 1
fi

READS_DIR=$3
if [ $READS_DIR != "" -a -d $READS_DIR ]
then
  echo "${READS_DIR} found :)"
else
  echo "${READS_DIR} not found :("
  exit 1
fi

STDOUT="${READS_DIR}/stdout.lisa"
echo $STDOUT

STDERR="${READS_DIR}/stderr.lisa"
echo $STDERR


for READS in `ls -v $READS_DIR/*.fq`
do
    echo "Processing reads file: $READS"
    RESULT="$READS.lisa"
    COMMAND="/usr/bin/time -p $LISA solve $DB $INDEX $READS $RESULT ${LISA_OPTS}"
    echo "Command: $COMMAND"
    $COMMAND
done > $STDOUT 2> $STDERR
