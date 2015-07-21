#!/bin/bash

rm -f temp.dat
rm -f sum.dat
rm -f col1.dat
touch temp.dat
touch sum.dat
touch col1.dat

for file in leb_orient*.dat
do 
awk  'BEGIN {ind=0}
     $0==" " {ind=ind+1; indNbr=NR}
     $1<30 && ind==0 {ypos[NR]=$2; print $1 >> "col1.dat"}
     $1<30 && ind==2 {yneg[NR-indNbr]=$2}
     $1==30 && ind==0 {linesNbr=NR-1}
     END {
     for(i=1; i<=linesNbr; i++) {
     print ypos[i] >> "temp.dat"
     } 
     printf("\n\n") >> "temp.dat"
     for(i=1; i<=linesNbr; i++) {
     print yneg[i] >> "temp.dat"
     } 
     print FILENAME
     } '   "$file"
if [ -s "sum.dat" ] 
then
  paste -d " " sum.dat temp.dat > _sum_.dat
  mv _sum_.dat sum.dat  
  rm temp.dat
else
  cp temp.dat sum.dat
  rm temp.dat
fi

done

touch temp.dat

awk 'BEGIN {sum=0}
     $0==" " {print("\n") >> "temp.dat"}
     $0!=" " {for(i=1; i<=NF; i++) sum=sum+$i; print sum >> "temp.dat"; sum=0;}' sum.dat
     
paste -d " " col1.dat temp.dat > _sum_.dat

mv _sum_.dat sum.dat  
rm temp.dat
rm col1.dat 

