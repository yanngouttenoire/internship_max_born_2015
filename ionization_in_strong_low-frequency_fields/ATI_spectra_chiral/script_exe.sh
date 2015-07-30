#!/bin/bash

for ((i=1; i<=38; i++))
do
qsub script_$i.pbs
sleep 1
done
