#!/bin/sh


awk '{printf "%f %f\n",$1, $2}' $1 > output

cp output $1
rm output 

