#!/bin/bash
# to set max min values on x y axis
# 


function xymaxmin(){

xmax=$(awk 'NR == 1{max = $1}{if($1 > max) max = $1}END{print max }' 1.txt) 
xmin=$(awk 'NR == 1{min = $1}{if($1 < min) min = $1}END{print min }' 1.txt) 
for i in $@
do
   xmaxn=$(awk 'NR == 1{max = $1}{if($1 > max) max = $1}END{print max }' $i) 
   xminn=$(awk 'NR == 1{min = $1}{if($1 < min) min = $1}END{print min }' $i) 
 

   if [ $(bc <<< "$xmax <= $xmaxn") -eq 1 ] ;then   xmax=$xmaxn   ;fi
   if [ $(bc <<< "$xminn <= $xmin") -eq 1 ] ;then   xmin=$xminn   ;fi

done

#difx=$(echo " scale=3;($xmax - $xmin)/10.0" | bc)
#xmax=$(echo " scale=3;$xmax + $difx" | bc)
#xmin=$(echo " scale=3;$xmin - $difx" | bc)




# for y axis ----------------------------

ymax=$(awk 'NR == 1{max = $2}{if($2 > max) max = $2}END{print max }' 1.txt) 
ymin=$(awk 'NR == 1{min = $2}{if($2 < min) min = $2}END{print min }' 1.txt)

for i in $@
do
   ymaxn=$(awk 'NR == 1{max = $2}{if($2 > max) max = $2}END{print max }' $i) 
   yminn=$(awk 'NR == 1{min = $2}{if($2 < min) min = $2}END{print min }' $i) 
 

   if [ $(bc <<< "$ymax <= $ymaxn") -eq 1 ] ;then   ymax=$ymaxn   ;fi
   if [ $(bc <<< "$yminn <= $ymin") -eq 1 ] ;then   ymin=$yminn   ;fi

done

dify=$(echo " scale=3;($ymax - $ymin)/10.0" | bc)
ymax=$(echo " scale=3;$ymax + $dify" | bc)
ymin=$(echo " scale=3;$ymin - $dify" | bc)



###-----main-----###

r=$xmin/$xmax/$ymin/$ymax
echo $r
} 

