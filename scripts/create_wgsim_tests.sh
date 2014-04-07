#za zadanu fasta bazu ova skripta generira single readove 
#razlicitih duljina i gresaka

# ovo je mozda ubuduce dobro parametrizirati ili 
# podrzati da uzima wgsim iz PATH-a
WGSIM=~/wgsim/wgsim
# WGSIM=~/samtools/misc/wgsim 
# WGSIM=/home/isovic/benchmark2/simulators/lh3-wgsim-a12da33/wgsim
                    

if [ $# -ne 2 ]
then
    echo "Usage: ./create_wgsim_tests.sh <path-to-fasta-db> <output-directory>"
    exit 1
fi

DB=$1
OUTPUT_DIR=$2
RANDOM_SEED=1603

if [ $DB != "" -a -f $DB ]
then
  echo "$DB found :)"
else
  echo "$DB not found :("
  exit 1
fi

mkdir -p $OUTPUT_DIR
echo "$OUTPUT_DIR found or created!"

for BASE_ERROR_RATE in 0.02 0.05 0.10
do
    for READ_LEN in 200 500 1000 2000 5000 10000
    do 
	echo "Creating test case with error=$BASE_ERROR_RATE and len=$READ_LEN"
	OUTPUT_FILE="$OUTPUT_DIR/error${BASE_ERROR_RATE}_len${READ_LEN}.fq"
	COMMAND="$WGSIM -N 100000 -e $BASE_ERROR_RATE -A 0.1 -d $READ_LEN -1 $READ_LEN $DB $OUTPUT_FILE /dev/null >/dev/null" # 2>/dev/null"
	echo "Command: $COMMAND"
	bash -c "$COMMAND"
    done
done
