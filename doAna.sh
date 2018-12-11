#!/bin/bash
case $1 in
  22) ANA=MC_TTBAR;;
  23) ANA=MC_TTBAR_SINGLEDECAY_23;;
  24) ANA=MC_TTBAR_SINGLEDECAY;;
  26) ANA=MC_TTBAR_DOUBLEDECAY;; 
esac
shift

outputs=()
for input in $@
do
  output="${input%.*}.yoda"
  if [[ $output == *"-"* ]]; then
    output="${output%-*}.yoda" 
  fi

  outputs+=($output)

  echo $output
  echo $input
  rivet --analysis $ANA ${input//\\} -o $output --quiet & 

done

