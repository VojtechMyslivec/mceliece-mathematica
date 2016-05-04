#!/bin/bash
typeset -i m t n r k;
for (( m=2 ; m<=8 ; m++ )); do
    printf "\n%2s %2s %3s %3s %3s\n" m t n r k

    for (( t=1 ; 1; t++ )); do 
        ((n=2**m))
        ((r=t*m))
        ((k=n-r))
        [ "$k" -gt 0 ] || break

        printf "%2d %2d %3d %3d %3d\n" "$m" "$t" "$n" "$r" "$k"
    done
done
