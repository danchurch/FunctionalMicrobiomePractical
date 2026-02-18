Read_set_list=(
"ERR3801502"
"ERR3801542"
"ERR3801603"
)

# Download reads from NCBI

for READ_SET in ${Read_set_list[@]};

do
  fasterq-dump \
    --progress \
    --split-files \
    --skip-technical \
    --outdir ${READ_SET} \
    ${READ_SET}
done


# Do a quick quality control of the reads (99% of reads are retained)

for READ_SET in ${Read_set_list[@]};

do
  bbduk.sh \
    threads=12 \
    ref=adapters,artifacts,phix \
    qtrim=rl \
    trimq=10 \
    maxns=0 \
    in=${READ_SET}/${READ_SET}_1.fastq \
    in2=${READ_SET}/${READ_SET}_2.fastq \
    out=${READ_SET}_1.fastq \
    out2=${READ_SET}_2.fastq \
    overwrite=T

done
