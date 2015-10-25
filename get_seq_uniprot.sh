#!/bin/bash

#download fasta seqs given file of uniprot ids

file=$1
name=$2

list=$(cat ${1})

#mkdir ${name}
cp ${2} ${name}
cd ${name}

for word in ${list}
do
    wget -nv http://www.uniprot.org/uniprot/$word.fasta
done
