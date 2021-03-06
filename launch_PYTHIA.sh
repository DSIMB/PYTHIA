#! /bin/bash

start=`date +%s`

# Show this for help
usage() {
cat << EOF
Usage:
        -i    | --seq         (Required)     Path to input Fasta sequence file. The path is relative to the project folder.
        -l    | --model       (Required)     Use "balanced" or "global" accuracy model.
        -d    | --database    (Required)     Name of the database for HHBlits.
                                             Ex: UniRef30_2020_06
                                             If your path is '/path/to/Uniclust/UniRef30_2020_03_a3m.ffdata'
					     If launching with Docker:
                                             	please provide: -d UniRef30_2020_03
					     If launching with Conda:
						please provide: -d /path/to/Uniclust/UniRef30_2020_03
        -o    | --outdir      (Required)     Path to the output directory.
        -c    | --cpus        (Optionnal)    Number of CPUs to use (for HHblits). Set to 0 for all available memory. Default is 0.
        -m    | --memory      (Optionnal)    Maximum RAM to use in Gb (for HHblits) Set to 0 for all. Default is 0.
        -h    | --help        (Optionnal)    Brings up this help
EOF
}


# Detect maximum hardware ressources based on OS
MAX_CPUS=$(getconf _NPROCESSORS_ONLN)
# get available memory (inactive RAM) in GB
AVAIL_MEMORY=$(awk '/MemFree/ { printf "%.0f \n", $2/1024/1024 }' /proc/meminfo)

# Set default resources based on HHblits defaults
CPUS=2
MEMORY=3


# Parse command line arguments
if [ "$#" -lt 4 ]; then
    printf "At least 4 parameters required.\n"
    usage
    exit 1
fi

while [ "$1" != "" ]; do
    case $1 in
        -i | --seq )
            shift
            SEQ=$1
        ;;
        -l | --model )
            shift
            MODEL=$1
        ;;
        -d | --database )
            shift
            DATABASE=$1
        ;;
        -o | --outdir )
            shift
            OUTDIR=$1
        ;;
        -c | --cpus )
            shift
            CPUS=$1
        ;;
        -m | --memory )
            shift
            MEMORY=$1
        ;;
        -h | --help )
            usage
            exit
        ;;
        * ) usage
            exit 1
    esac
    shift
done


# Check if we are living in Docker container or not
if grep -q docker /proc/1/cgroup; then
    PROJECT_DIR=/project
    DATABASE_DIR=/database
else
    PROJECT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
    DATABASE_DIR=$(dirname -- $DATABASE)
    DATABASE=$(basename -- $DATABASE)
fi


# Check for valid sequence input
if [ ! -f $PROJECT_DIR/$SEQ ]; then
    printf "\nA valid input sequence file is required, provide it with the flag: -i sequence.fasta\nMake sure the path to the file is relative to your project folder.\n\n"
    exit
elif [ $(grep -c "^>" $PROJECT_DIR/$SEQ) -gt 1 ]; then
    printf "\nPlease provide only one FASTA sequence.\n\n"
    exit
elif [ $(grep -c "^>" $PROJECT_DIR/$SEQ) -eq 0 ]; then
    printf "\nPlease provide at least one FASTA sequence.\n\n"
    exit
fi


TMP_SEQ=/tmp/seq.fasta
# Convert multiline fasta to single line
sed ':a;N;/^>/M!s/\n//;ta;P;D' < $PROJECT_DIR/$SEQ > $TMP_SEQ
# Extract only sequence without header (and trim whitespaces on extremities)
SEQ=$(sed ':a;N;/^>/M!s/\n//;ta;P;D' < $PROJECT_DIR/$SEQ | tail -1 | xargs)
printf "\nInput sequence:\n$SEQ\n\n"
# Check if sequence contains whitespaces
if [[ $SEQ =~ .*[[:space:]]+.* ]]; then
    printf "\nThe sequence contains whitespaces, please provide a valid sequence.\n\n"
    exit
# Check unsupported amino acids
elif [[ $SEQ =~ [UXBZ]+ ]]; then
    printf "\nYour sequence contains unsupported amino acids (U, X, B or Z).\n\n"
    exit
# Check if sequence if a correct protein sequence
elif [[ ! $SEQ =~ ^[ACDEFGHIKLMNPQRSTVWY]+$ ]]; then
    printf "\nPlease provide a valid input protein FASTA sequence file.\n\n"
    exit
# Check if it is a DNA sequence
elif [[ $SEQ =~ ^[atgcATGC]+$ ]]; then
    printf "\nYour sequence looks like a DNA sequence.\nPlease provide a protein sequence.\n\n"
    exit
fi

# Check if database is accessible
if [ ! -f $DATABASE_DIR/$DATABASE"_a3m.ffindex" ]; then
    printf "\nThe database $DATABASE_DIR/$DATABASE seem unusable, please link the right database."
    printf "\nIf your path is '/path/to/database/UniRef30_2020_03_a3m.ffdata' provide: -d UniRef30_2020_03 when using Docker, or -d /path/to/database/UniRef30_2020_03 when using Conda.\n\n"
    exit
fi

# Check if the output directory exists
if [ ! -d $PROJECT_DIR/$OUTDIR ]; then
    printf "\n$OUTDIR does not exist, please provide a right output directory.\n\n"
    exit
fi

# Check for valid ressources args (integers)
if [[ ! $CPUS =~ ^[0-9]+$ || $CPUS -gt $MAX_CPUS ]]; then
    printf "\nThe number of CPUs argument should be an integer 1 >= cpus >= $MAX_CPUS.\n\n"
    exit
elif [[ $CPUS -eq 0 ]]; then
    CPUS=$MAX_CPUS
fi

if [[ ! $MEMORY =~ ^[0-9]+$ ]]; then
    printf "\nThe memory argument should be an integer 3 >= memory >= $AVAIL_MEMORY.\n\n"
    exit
elif [[ $MEMORY -eq 0 ]]; then
    MEMORY=$AVAIL_MEMORY
elif [[ $AVAIL_MEMORY -lt 3 ]]; then
    printf "\nOnly $AVAIL_MEMORY Gb of RAM is available.\nHHblits recommands at least 3 GB of memory, the program might not run well."
    printf "\nSetting minimum 1 Gb of memory.\n"
    if [[ "$OSTYPE" == "darwin"* ]]; then
        printf "\nPlease consider increasing the runtime memory allocated to the Docker Desktop program on your Mac.\n"
    elif [[ "$OSTYPE" == "cygwin" ]]; then
        printf "\nPlease consider increasing the runtime memory allocated to the Docker Desktop program on your Windows PC.\n"
    fi
    printf "\n"
elif [[ $MEMORY -lt 3 && $AVAIL_MEMORY -gt 3 ]]; then
    printf "\nYou asked for $MEMORY Gb of RAM."
    printf "\nHHblits recommands at least 3 GB of memory, the program might not run well."
    printf "\nYou have $AVAIL_MEMORY GB available so we will use default 3 Gb instead.\n\n"
    MEMORY=3
elif [[ $MEMORY -lt 3 && $AVAIL_MEMORY -lt 3 && $MEMORY -lt $AVAIL_MEMORY ]]; then
    printf "\nYou asked for $MEMORY Gb of RAM."
    printf "\nYou have $AVAIL_MEMORY GB available."
    printf "\nSetting minimum value of 1 Gb instead."
    printf "\nHHblits recommands at least 3 GB of memory, the program might not run well.\n\n"
    MEMORY=1
elif [[ $MEMORY -gt $AVAIL_MEMORY ]]; then
    printf "\nOnly $AVAIL_MEMORY GB memory available, you asked for $MEMORY GB."
    printf "\nUsing $AVAIL_MEMORY GB memory instead.\n\n"
    MEMORY=$AVAIL_MEMORY
fi

printf "Using:\n  $MEMORY Gb max. of RAM\n  $CPUS cpus\n\n"


# Check correct type of model asked
if [[ $MODEL == "global" ]]
then
    MODEL_ARCH="models/global_accuracy/global_acc_model_arch.json"
    MODEL_WEIGHTS="models/global_accuracy/global_acc_model_weights.hdf5"
elif [[ $MODEL == "balanced" ]]
then
    MODEL_ARCH="models/balanced_accuracy/balanced_acc_model_arch.json"
    MODEL_WEIGHTS="models/balanced_accuracy/balanced_acc_model_weights.hdf5"
fi

if [[ ! -f $MODEL_ARCH ]]
then
    printf "\nUnable to find $MODEL_ARCH\n"
    exit 1
elif [[ ! -f $MODEL_WEIGHTS ]]
then
    printf "\nUnable to find $MODEL_WEIGHTS\n"
    exit 1
fi

# Create a directory for the job output
JOBDIR=$(mktemp -d -t PYTHIA.XXXX --suffix=-$(date +%Y%m%d%H%M%S) -p $PROJECT_DIR/$OUTDIR)
ORIGINAL_OUTDIR_NAME=$OUTDIR/`basename $JOBDIR`

# Set global paths

# The project is mounted as docker bind volume
# into $PROJECT_DIR directory in the docker container
PROJECT=$PROJECT_DIR
HHSUITE=$PROJECT/bin/hh-suite
HHBLITS=$HHSUITE/bin/hhblits
HHFILTER=$HHSUITE/bin/hhfilter
# path to database is a mounted volume to /database in Docker image env
# the name of the database is given in command line argument
DBHHBLITS=$DATABASE_DIR/$DATABASE
OUTDIR=$JOBDIR
SEQ=$TMP_SEQ
SCRIPTS=$PROJECT/src
DATA=$PROJECT/data
MODELS=$PROJECT/models
WINDOW=61



#############
### START ###
#############


### Step 1. Create pssm data.
### HHblits is time consuming, probably remove the -realign_old_hits option. Probably decrease -max_filt and -realign_max.
printf "Run HHblits ... "
$HHBLITS -cpu $CPUS -maxmem $MEMORY -maxfilt 10000 -diff inf -B 10000 -Z 10000 -e 0.0001 -cov 75 -realign_old_hits -realign_max 10000 -n 3 -i $SEQ -d $DBHHBLITS -oa3m $OUTDIR/job.a3m -o $OUTDIR/job.hhr 1>/dev/null 2> $OUTDIR/hhblits.log
printf "done\n"

printf "Run HHfilter ... "
$HHFILTER -v 2 -id 99 -neff 20 -qsc -30 -cov 75 -i $OUTDIR/job.a3m -o $OUTDIR/job\_filter.a3m 2>/dev/null
printf "done\n"

#perl $HHSUITE/scripts/reformat.pl -M -uc first a3m fas job\_filter.a3m pdb_mfasta/job.mfasta
# without insertions in the first (query) sequence
printf "Reformat ... "
perl $HHSUITE/scripts/reformat.pl -r -M -uc first a3m fas $OUTDIR/job\_filter.a3m $OUTDIR/job.mfasta_woi 1>/dev/null 2>&1
printf "done\n"

printf "Transform alignment to frequence matrix ... "
$SCRIPTS/ali2freq-py3.py -first -al $OUTDIR/job.mfasta_woi -m $DATA/homstradfreq.txt -gapaa 1> $OUTDIR/job.aamtx_gaps 2>/dev/null
printf "done\n"

printf "Create AAindex and one-hot encodings ... "
### Step 2. Create aaindex data.
$SCRIPTS/fasta2vector_wgap.pl $SEQ $DATA/selected_aaindex1_reformated_58_Van_Westen > $OUTDIR/job.aaindex

### Step 3. Create one-hot data.
$SCRIPTS/seq2onehot.pl $SEQ AA >> $OUTDIR/job.onehot
printf "done\n"

### Step 4. Translate each encoding to a vector of the given window size.
printf "Translate encodings into vectors ... "
$SCRIPTS/prepare_data_for_learning.pl $OUTDIR/job.aamtx_gaps $WINDOW > $OUTDIR/job.vector_aamtx
$SCRIPTS/prepare_data_for_learning.pl $OUTDIR/job.aaindex $WINDOW > $OUTDIR/job.vector_aaindex
$SCRIPTS/prepare_data_for_learning.pl $OUTDIR/job.onehot $WINDOW > $OUTDIR/job.vector_onehot
printf "done\n"

### Step 5. Merge all features to one vector
printf "Merge vectors ... "
python3 $SCRIPTS/create_vector_features.py -i $OUTDIR/job -o $OUTDIR/job.merged -w $WINDOW
printf "done\n"


### Step 6. Make prediction
# Activate Tensorflow conda environment to run the script
printf "Run prediction ... "
$SCRIPTS/pythia.py -m $OUTDIR/job.merged -o $OUTDIR/prediction -f $SEQ -mjson $MODEL_ARCH -mh5 $MODEL_WEIGHTS 1>$OUTDIR/predictions.log 2>&1
printf "done\n\n"

# Create downloadable tar.gz
mkdir $OUTDIR/logs
cp $OUTDIR/job* $OUTDIR/*.log $OUTDIR/logs/
tar -czf $OUTDIR/pythia_job\_results.tar.gz $OUTDIR/prediction/*.csv $OUTDIR/prediction/*.fasta $OUTDIR/logs 1>/dev/null 2>&1
rm $OUTDIR/job* $OUTDIR/*.log
printf "Results can be found in $ORIGINAL_OUTDIR_NAME\n\n"
end=`date +%s`
runtime=$((end-start))
printf "Total runtime: $runtime seconds\n\n"

