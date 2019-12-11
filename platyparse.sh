#!/bin/bash

perl -lane '
  BEGIN{ print join("\t","CHROM","POS","REF","ALT","FILTER","TC","TR","FREQREF","FREQALT","MAF") }
  if(length($F[3])==1 && length($F[4])==1){
      $TC = $1 if $F[7] =~ /TC\=(\d+)/;
      $TR = $1 if $F[7] =~ /TR\=(\d+)/;
      print join("\t",$F[0],$F[1],$F[3],$F[4],$F[6],$TC,$TR,(($TC-$TR)/$TC),($TR/$TC),(($TC-$TR)/$TC)<($TR/$TC) ? (($TC-$TR)/$TC) : ($TR/$TC));
  }' $1
