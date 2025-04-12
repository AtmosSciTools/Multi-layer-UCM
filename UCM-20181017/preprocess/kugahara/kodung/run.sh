

fname="kugahara_0901_"
ftail=".txt"

for var in sd lw temp qv u ah 
do 
name=$fname$var$ftail
echo $name

awk '{print $2}' $name > $var

done

paste sd lw temp qv u ah > $fname"vars"$ftail

rm sd lw temp qv u ah


