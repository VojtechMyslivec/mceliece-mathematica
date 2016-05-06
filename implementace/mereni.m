(* ::Package:: *)

(* ::Text:: *)
(*Vojt\[EHacek]ch Myslivec, FIT \[CapitalCHacek]VUT v Praze, 2015/2016 *)


(* ::Title:: *)
(*Mereni*)


(* ::Section:: *)
(*Casov\[EAcute] slozitosti*)


(* ::Subsection:: *)
(*Zavislosti*)


SetDirectory[NotebookDirectory[]];
Get["src/moduly.m"]
Get["src/rozsirenaBinarniTelesa.m"]
Get["src/ireducibilniBinarniGoppaKody.m"]
Get["src/mceliece.m"]


(* ::Subsection:: *)
(*Mereni*)


pocetMereni = 10;
mereniM = {};

For[ m = 3, m <= 8, m++,
    mereniT = {};

    For[ t = 2, 2^m - m*t > 0, t++,
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
        soubor = StringTemplate[ "mereni/data_m`1`_t`2`.txt" ][m,t];
        Export[ soubor, {list} ];

        AppendTo[ mereniT, list ];
    ];
    AppendTo[ mereniM, { m, mereniT } ];
]


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


(* ::Subsection:: *)
(*Vykresleni grafu*)


(* ::Text:: *)
(*Ziskani dat pro grafy*)


mereniT = mereniM[[3,2]];
Transpose[ { mereniT[[All,1]], mereniT[[All,2,1,1]] } ] // MatrixForm
Transpose[ { mereniT[[All,1]], mereniT[[All,2,2,1]] } ] // MatrixForm
Transpose[ { mereniT[[All,1]], mereniT[[All,2,3,1]] } ] // MatrixForm


(* ::Text:: *)
(*Grafy*)


mereniT = mereniM[[6,2]];
listGen    = Transpose[ { mereniT[[All,1]], mereniT[[All,2,1,1]] } ];
listSifr   = Transpose[ { mereniT[[All,1]], mereniT[[All,2,2,1]] } ];
listDesifr = Transpose[ { mereniT[[All,1]], mereniT[[All,2,3,1]] } ];

ListPlot[ listGen,    AxesOrigin->0 ]
ListPlot[ listSifr,   AxesOrigin->0 ]
ListPlot[ listDesifr, AxesOrigin->0 ]


(* ::Text:: *)
(*Tabulka hodnot: { t, generovani, sifrovani, desifrovani }*)


Transpose[ { mereniT[[All,1]], mereniT[[All,2,1,1]], mereniT[[All,2,2,1]], mereniT[[All,2,3,1]] } ] // MatrixForm
