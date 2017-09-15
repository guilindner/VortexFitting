for f in *.pdf
do
  echo "convert $f ${f%%.*}.png"
  convert $f ${f%%.*}.png
done
