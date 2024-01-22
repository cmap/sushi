DOCKER_TAG="develop"

#### fastq2readcounts ###
#docker run -it -v $PWD/assets/:/data prismcmap/sushi:$DOCKER_TAG fastq2readcount \
#    --fastq /data/fastq/ \
#    --out /data/out \
#    -i1 _I1_ -i2 _I2_ -b _R1_

FOUND_LINE=$(cat ./assets/out/raw_counts.csv | grep -c "ACAGGATG,AAGTAGAG,ACATTACTTCCATATACAACTAAT,49")

if [ "$FOUND_LINE" -eq "1" ]; then
    echo "fastq2readcounts test passed"
    rm -f ./assets/out/
else
    echo "fastq2readcounts test failed"
    exit 1
fi

