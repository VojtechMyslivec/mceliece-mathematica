#!/bin/bash
# spoji data z mereni.m v adresari $1

set -o pipefail
set -o nounset

USAGE="USAGE
    $0 adresar_s_daty

    V adresari s daty ocekava data z mereni z notebooku mereni.m
    se jmeny ve formaty data_m{m}_t{t}.txt

    Data z techto souboru spoji a vypise ve formatu pro software
    mathematica."


[ $# -eq 1 ] || {
    echo "$USAGE" >&2
    exit 1
}

[ "$1" == "-h" -o "$1" == "--help" -o "$1" == "help" ] && {
    echo "$USAGE"
    exit 0
}

dir=$1

[ -d "$dir" -a -r "$dir" ] || {
    echo "Argument '$dir' neni citelny adresar." >&2
    exit 1
}

# pattern pro jmeno souboru, pozor na RE
maska_prefix="data_m"
maska_stred="_t"
maska_suffix=".txt"

soubory=$( ls "$dir/$maska_prefix"*"$maska_stred"*"$maska_suffix" 2>&1 ) || {
    echo "V adresari '$dir' nejsou data v ocekavanem souboru." >&2
    exit 2
}

# serazene hodnoty m
hodnotyM=$(
    sed "s|^${dir}/${maska_prefix}||; s|${maska_stred}.*${maska_suffix}$||" \
        <<< "$soubory" \
        | sort -nu
)
posledniM=$(
    tail -1 <<< "$hodnotyM"
)

# zacatek vystupu -- seznamu
echo "mereniM = {"
while read m; do

    # serazene hodnoty t pro dane m
    hodnotyT=$(
        sed -n "s|${maska_suffix}$||
                s|^${dir}/${maska_prefix}${m}${maska_stred}||p" \
            <<< "$soubory" \
            | sort -nu
    )
    posledniT=$(
        tail -1 <<< "$hodnotyT"
    )

    # debugovaci vypis
    echo "$m:" $hodnotyT >&2

    # zacatek prvku pro m
    echo    "  {"
    echo    "    ${m},"
    echo    "    {"

    # hodnoty t
    while read t; do
        # soubory neobsahuji newline, proto trapeni echa
        echo -n "$( cat "$dir/data_m${m}_t${t}.txt" | sed 's|^|      |' )"
        if [ "$t" != "$posledniT" ]; then
            echo ","
        else
            echo
        fi
    done <<< "$hodnotyT"

    # konec prvku pro m
    echo    "    }"
    echo -n "  }"
    # osetreni chybejici carky za poslednim prvkem
    if [ "$m" != "$posledniM" ]; then
        echo ","
    else
        echo
    fi

done <<< "$hodnotyM"

# konec vystupu -- seznamu
echo "};"

