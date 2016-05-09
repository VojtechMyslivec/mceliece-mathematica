(* ::Package:: *)

(* ::Text:: *)
(*Vojtech Myslivec, FIT CVUT v Praze, 2015/2016 *)


(* ::Title:: *)
(*Mereni*)


(* ::Section:: *)
(*Mereni*)


(* ::Subsection:: *)
(*Zavislosti*)


SetDirectory[NotebookDirectory[]];
Get["src/moduly.m"]
Get["src/rozsirenaBinarniTelesa.m"]
Get["src/ireducibilniBinarniGoppaKody.m"]
Get["src/mceliece.m"]


(* ::Subsection:: *)
(*Mereni*)


(* ::Text:: *)
(*Hodnoty jsou ve formatu*)
(*{*)
(* { *)
(*   m,*)
(*   {*)
(*     { *)
(*       t,*)
(*       {*)
(*         { tGen, { tMatV, tMatD, tMatH, tMatG, tSyndromy, tMatS, tMatP, tMatHatG, tInv } }*)
(*         { tSifr }*)
(*         { tDesifr, { tPolyS, tPolyR, tEEA, tPolySigma, tDosazeni } }*)
(*       }*)
(*     },*)
(*     ...*)
(*   }*)
(* },*)
(* ...*)
(*}*)


mMin = 3; mMax = 8; tMin = 2;


pocetMereni = 10;
mereniM = {};

For[ m = mMin, m <= mMax, m++,
    mereniT = {};

    For[ t = tMin, 2^m - m*t > 0, t++,
        Print[ "m = ", m, ", t = ", t ];

        sumaGen = sumaSifr = sumaDesifr = 0;
        pocetUspesnychMereni = 0;
        For[ i = 0, i < pocetMereni, i++,

            (* Generovani kodu *)
            (* tKod : { tMatV, tMatD, tMatH, tMatG, tSyndromy } *)
            (* tMc  : { tMatS, tMatP, tMatHatG, tInv } *)
            tGen = AbsoluteTiming[
                genKod = Quiet[ generujMcEliece[ m, t]];
            ][[1]];
            If[ genKod == Null,
                Continue[];
            ];
            pocetUspesnychMereni++;

            { { tKod, tMc }, { soukromyKlic, verejnyKlic, parametry } } = genKod;
            sumaGen += { tGen, Flatten[ { tKod, tMc } ] };

            (* Sifrovani *)
            zprava = nahodnyPolynom[ 2, parametry[[2]] ];
            tSifr = AbsoluteTiming[
                c = sifrujMcEliece[ zprava, verejnyKlic, parametry ];
            ][[1]];
            sumaSifr += {tSifr};

            (* Desifrovani *)
            (* tDek : { tPolyS, tPolyR, tEEA, tPolySigma, tDosazeni } *)
            tDesifr = AbsoluteTiming[
                { tDek, mm } = desifrujMcEliece[ c, soukromyKlic, parametry ];
            ][[1]];

            If[ mm != zprava,
                Print[ "Stalo se neco spatneho!?" ];
            ];

            sumaDesifr += { tDesifr, tDek };
        ];

        If[ pocetUspesnychMereni == 0,
            (* then *)
            Print[ "Nepodarilo se namerit" ];
            Continue[];,

            (* else *)
            Print[ "Pocet uspesnych mereni: ", pocetUspesnychMereni ];
        ];

        { prumerGen, prumerSifr, prumerDesifr } = { sumaGen, sumaSifr, sumaDesifr }/pocetUspesnychMereni;

        (* export dat do souboru *)
        list = { t, { prumerGen, prumerSifr, prumerDesifr } };
        soubor = StringTemplate[ "data/data_m`1`_t`2`.txt" ][m,t];
        Export[ soubor, {list} ];

        AppendTo[ mereniT, list ];
    ];
    AppendTo[ mereniM, { m, mereniT } ];
]

