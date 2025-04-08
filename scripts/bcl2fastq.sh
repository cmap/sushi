#!/bin/bash
source /broad/software/scripts/useuse
reuse .bcl2fastq2-v2.20.0 > /dev/null
reuse .bcl2fastq2-2.20.0.422 > /dev/null

#OUT_DIR=/xchip/prism/bcl2fastq/
PSEQ_obelix=/cmap/obelix/pod/prismSeq/
#OUT_bcl2fastq=/xchip/prism/bcl2fastq/ #new

#optional
if test $# -lt 1; then
  printf "Usage ./bcl2fastq.sh [options]\nOptions include:\n"
  printf -- "-s, --seq_code \t\t Sequencer run code E.g. JNV9V \n"
  printf -- "-p, --proj_code \t Project code to prep project directory in /cmap/obelix/pod/prismSeq/ \n"
  printf -- "-b, --build_dir \t Build directory, usually on /cmap/obelix/. Overrides PROJ_CODE\n"
  printf -- "-o, --out_dir \t\t Path to temp storage of fastq files on /xchip/prism/ \n"
  printf -- "-r, --runfolder_dir \t Path to folder with BCL files \n" #new
  printf -- "-h, --help \t\t Print this help text\n"
  exit 1
fi

while test $# -gt 0; do
  case "$1" in
    -h|--help)
      printf "Usage ./bcl2fastq.sh [options]\nOptions include:\n"
      printf -- "-s, --seq_code \t\t Sequencer run code E.g. JNV9V \n"
      printf -- "-p, --proj_code \t Project code to prep project directory in /cmap/obelix/pod/prismSeq/ \n"
      printf -- "-b, --build_dir \t Build directory, usually on /cmap/obelix/. Overrides PROJ_CODE\n"
      printf -- "-o, --out_dir \t\t Path to temp storage of fastq files on /xchip/prism/ \n"
      printf -- "-r, --runfolder_dir \t Path to folder with BCL files \n" #new
      printf -- "-h, --help \t\t Print this help text\n"
      exit 0
      ;;
    -s|--seq_code)
      shift
      SEQ_CODE=$1
      ;;
    -b|--build_dir)
      shift
      #echo $1
      BUILD_DIR=$1
      ;;
    -p|--proj_code)
      shift
      #echo $1
      PROJ_CODE=$1
      ;;
    -o|--out_dir)
      shift
      #echo $1
      OUT_DIR=$1
      ;;
    #new
    -r|--runfolder_dir)
      shift
      #echo $1
      RUNFOLDER_DIR=$1
      ;;
    *)
      printf "Unknown parameter: %s \n" "$1"
      shift
      ;;
  esac
  shift
done

#new
if [ -z $OUT_DIR ]
then
  OUT_DIR=$(echo /xchip/prism/bcl2fastq/$SEQ_CODE)
fi

if [ ! -d $OUT_DIR ]
then
  mkdir $OUT_DIR
fi

#new
if [ -z $RUNFOLDER_DIR ]
then
  RUNFOLDER_DIR=$(echo /xchip/prism/MiSeq\ Outputs/*-$SEQ_CODE)
fi

#RUNFOLDER_DIR=$(echo /xchip/prism/MiSeq\ Outputs/*-$SEQ_CODE)

echo $RUNFOLDER_DIR

#bcl2fastq --runfolder-dir "$RUNFOLDER_DIR" --output-dir $OUT_DIR/$SEQ_CODE --minimum-trimmed-read-length 35 --mask-short-adapter-reads 22 --create-fastq-for-index-reads
bcl2fastq --runfolder-dir "$RUNFOLDER_DIR" --output-dir $OUT_DIR --minimum-trimmed-read-length 35 --mask-short-adapter-reads 22 --create-fastq-for-index-reads

if [ -z $BUILD_DIR ]
then
  BUILD_DIR=$PSEQ_obelix/$PROJ_CODE
fi

if [ ! -d $BUILD_DIR ]
then
  mkdir $BUILD_DIR
fi

if [ ! -d $BUILD_DIR/fastq/ ]
then
  mkdir -p -m a=rwx $BUILD_DIR/fastq/
fi

#echo Copying fastq files from $OUT_DIR/$SEQ_CODE/ to $BUILD_DIR/fastq/
#cp $OUT_DIR/$SEQ_CODE/*.fastq.gz $BUILD_DIR/fastq/

echo Copying fastq files from $OUT_DIR/ to $BUILD_DIR/fastq/ #new
cp $OUT_DIR/*.fastq.gz $BUILD_DIR/fastq/ #new
