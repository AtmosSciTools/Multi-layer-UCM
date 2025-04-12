#!/bin/sh


awk '{print NR, $1}' ah.txt


awk '$0 !~/^$/{print $2}' data_072930.txt > v
awk '$0 !~/^$/{print $3}' data_072930.txt > t

awk '$0 !~/^$/{print $2*1000000/3600}' data_sw072930.txt > sw

awk 'BEGIN{for(i = 0; i < 97; i++) print "-99"}' > hum
awk 'BEGIN{for(i = 0; i < 97; i++) print "-99"}' > lw
awk 'BEGIN{for(i = 0; i < 97; i++) print "-99"}' > sd






