names=(
"ERR3801502"
"ERR3801542"
"ERR3801603"
)

for i in ${names[@]}; do
  prefetch $i
  fasterq-dump --split-files $i
done
