#!/bin/bash

DIR=lebedev_38

rm -rf $DIR
mkdir -p $DIR
cp lebedev_table $DIR/.
cp script_exe.sh $DIR/.
make clean
make mrproper

for ((i=1; i<=38; i++))
do
sed "s/__myOrientation__/$i/" main.cpp.GEN > main.cpp
make
mv exec.out $DIR/exec_$i.out
sed "s/__EXE__/exec_$i.out/" script.pbs.GEN > $DIR/script_$i.pbs
echo "File n° $i"
done
