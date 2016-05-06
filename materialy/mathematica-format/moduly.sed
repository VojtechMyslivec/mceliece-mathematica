#!/bin/sed -f
# alespon trochu zformatuj mathematica package, ktery je ulozeny z nb formatu
# a tedy absolutne bez mezer a naplacato...

# vse e provadi nad nezakomentovanymi radky
# mezera kolem operatoru = <= >= ...
/^[ \t]*(\*/!s|\([^ ]\)\([<>!=:]\{0,2\}=\)\([^ ]\)|\1 \2 \3|g
# problemovy vyskyt
/^[ \t]*(\*/!s|]:=|] :=|g

#mezera kolem : ne ale :: a := a :>
/^[ \t]*(\*/!s|\([^: ]\):\([^>:= ]\)|\1 : \2|g


# moduly vyjimka!
# mezera za carkou a slozenou zavorkou ',' -> ', '
/^[ \t]*(\*/!s|\([{,]\)\([^ 0-9]\)|\1 \2|g
# mezera pred zaviraci slozenou zavorkou
/^[ \t]*(\*/!s|\([^0-9 ]\)}|\1 }|g

# mezera za hranatou zavorkou, krome [[ a \[
/^[ \t]*(\*/!s|\([^\[]\)\[\([^[ ]\)|\1[ \2|g
# a ]]
/^[ \t]*(\*/!s|\([^] ]\)\]\([^]]\)|\1 ]\2|g
# oprava \[Alpha ] apod.
/^[ \t]*(\*/!s|\(\\\[[^]]\+\) ]|\1]|g


# mezery kolem /@ nebo /. nebo /; nebo // nebo :>
/^[ \t]*(\*/!s,\([^ ]\)\(/@\|/\.\|/;\|//\|:>\)\([^ ]\),\1 \2 \3,g

# moduly vyjimky
/^[ \t]*(\*/!s|{i,2}|{ i, 2 }|g
/^[ \t]*(\*/!s|{2,\([0-9]\+\)},\([0-9]\+\)|{ 2, \1 }, \2|g

