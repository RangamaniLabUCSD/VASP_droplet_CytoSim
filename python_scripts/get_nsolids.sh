#!/bin/bash
while getopts f: flag
do
    case "${flag}" in
        f) filepath=${OPTARG};;
        #a) age=${OPTARG};;
        #f) fullname=${OPTARG};;
    esac
done
filename="${filepath}/objects.cmo"
startl=$(grep -n -a "section solid" $filename|head -1 | cut -d : -f 1)
endlA=$(grep -n -a "section single F" $filename|head -1 | cut -d : -f 1)
endlB=$(grep -n -a "section single A" $filename|head -1 | cut -d : -f 1)
endl=$endlA
if [ $endlA -gt $endlB] && [$endlB -gt 0]
   then endl=$endlB
fi
nline=$(grep -n -a "d1:" $filename|head -1| cut -d ' ' -f 2)
var="$((nline + 1))"
nsolids="$(((endl-startl-1)/var))"
#echo $startl
#echo $endlA
#echo $endlB
#echo $endl
#echo $nline
echo $nsolids
