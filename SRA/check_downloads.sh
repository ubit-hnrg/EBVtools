#!/bin/bash
runidfile=$1
cat $runidfile|while read s;
do
  echo $s;expected=0;obtained=1;
  expected=$(grep $s ids20200318.tsv |cut -f 16);
  obtained=$(zgrep '^@ER' $s/*_1*gz | wc -l);
  if [ $expected = $obtained ] ;
    then 
    echo $s'|'$expected'|'$obtained >> checked_ok.txt;
  else
  echo $s'|'$expected'|'$obtained >> doNotPassCheck.txt;
  fi;
done
