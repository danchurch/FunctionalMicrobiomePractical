names=(
"ERR3801502"
"ERR3801542"
"ERR3801603"
)

for READ_SETS in ${names[@]}; do
  fasterq-dump \
    --progress \
    --split-files \
    --skip-technical \
    --outdir ${READ_SETS} \
    ${READ_SETS}
done
