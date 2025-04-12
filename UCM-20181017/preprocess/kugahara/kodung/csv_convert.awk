#!/bin/bash

 awk '{out=$1; for(i=2;i<=NF;i++){out=out","$i}; print out}' $1
 
