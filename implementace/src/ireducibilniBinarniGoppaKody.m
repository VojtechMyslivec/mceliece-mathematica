(* ::Package:: *)

(* ::Text:: *)
(*Vojtech Myslivec, FIT CVUT v Praze, 2015/2016 *)


(* ::Title:: *)
(*Ireducibilni binarni Goppa kody*)


(* ::Subsection:: *)
(*Zavislosti*)


(* ::Text:: *)
(*Nutno pridat nasledujici "baliky"*)
(*"src/rozsirenaBinarniTelesa.m"*)
(*"src/moduly.m"*)


(* ::Section:: *)
(*GoppaKod -- implementace*)


(* ::Subsection::Closed:: *)
(*dotNad2*)


(* ::Text:: *)
(*Maticove nasobeni modulo 2*)


Unprotect[ dotNad2 ];
ClearAll[ dotNad2 ];

dotNad2[ A__ ] :=
    Mod[ Dot[A], 2 ]

Protect[ dotNad2 ];


(* ::Subsection::Closed:: *)
(*dotNadFF*)


(* ::Text:: *)
(*Maticove nasobeni nad telesem. TODO -- implementovat pomoci Inner[]?*)


Unprotect[ dotNadFF ];
ClearAll[ dotNadFF ];

dotNadFF::nkomp = "Nekompatibilni matice";

dotNadFF[
        A_List /; ArrayDepth[A] == 3,
        B_List /; ArrayDepth[B] == 3,
        M_List /; Length[M] == 2
    ] :=
    Module[ {
        da = Take[ Dimensions[A], 2 ], db = Take[ Dimensions[B], 2 ],
        matice = {}, radek, prvek, nula,
        i, j, k, p = NestWhile[ Last, M, ListQ ]
    },
    If[ da[[2]] != db[[1]],
        Message[ dotNadFF::nkomp ];
        Return[]
    ];
    nula = nulovyPolynom[ p, Dimensions[ A[[1,1]]] ];
    For[ i = 1, i <= da[[1]], i++,
        radek = { };
        For[ j = 1, j <= db[[2]], j++,
            prvek = nula;
            For[ k = 1, k <= da[[2]], k++,
                prvek = plus[ krat[ A[[i,k]], B[[k,j]], M ], prvek, p]
            ];
        AppendTo[ radek, prvek]
        ];
    AppendTo[ matice, radek]
    ];
    matice
]

Protect[ dotNadFF ];


(* ::Subsection::Closed:: *)
(*N-tice L -- pro ireducibilni Goppa kody*)


(* ::Text:: *)
(*Generuje n-tici L -- support Goppa kodu. Jsou to (permutovane) vsechny prvky telesa GF(2^m) a zadny z nich nesmi byt korenem polynomu g z GF(2^m)[x] (g je tedy ireducibilni).*)


Unprotect[ generujNticiL ];
ClearAll[ generujNticiL ];

(* vsechny kombinace, nahodne zamichane *)
generujNticiL[
        m_Integer
    ] :=
    RandomSample[ Tuples[ { 0, 1 }, m ] ];

Protect[ generujNticiL ];


(* ::Subsection::Closed:: *)
(*Matice D*)


(* ::Text:: *)
(*Matice D je soucast kontrolni matice Goppa kodu. Je to diagonalni matice nad GF(2^m), kde prvky na diagonale jsou 1/g(Subscript[L, i]). *)
(**)
(*     / 1/g(L1) 0       0       ... 0       \ *)
(*     | 0       1/g(L2) 0       ... 0       | *)
(*     | 0       0       1/g(L3) ... 0       | *)
(* D = | .       .       .       .   .       | *)
(*     | .       .       .        .  .       | *)
(*     | .       .       .         . .       | *)
(*     \ 0       0       0       ... 1/g(Ln) / *)
(**)


Unprotect[ maticeD ];
ClearAll[ maticeD ];

maticeD::ired = "Polynom g neni ireducibilni.";

maticeD[
        L_List,
        M : { g_List, MVnitrni : { ired_List, p_Integer }}
    ] :=
    Module[ {
        ll = Length[L],
        dosazeni, pravidla, nula
    },
    nula = nulovyPolynom[ p, Dimensions[ L[[1]]]];

    (* inverze prvku L dosazenych do polynomu *)
    (* TODO: vypocet dosadDoPolynomu lze zrychlit pomoci algoritmu Chienovo hledani korenu *)
    dosazeni = Quiet[ inverze[ dosadDoPolynomu[ g, #, MVnitrni ], MVnitrni ]& /@ L ];
    (* pokud nektera inverze neexistuje, neni to ireducibilni polynom *)
    If[ Length[ Cases[ dosazeni, Null ] ] != 0,
        Message[ maticeD::ired ];
        Return[]
    ];

    pravidla = Transpose[{ Range[ll], dosazeni }];
    pravidla = #[[1]]->#[[2]]& /@ pravidla;

    DiagonalMatrix[ Range[ll] ] /. {0->nula} /. pravidla
]

Protect[ maticeD ];


(* ::Subsection::Closed:: *)
(*Matice V*)


(* ::Text:: *)
(*Matice V je soucast kontrolni matice Goppa kodu. Je to Vandermondova matice nad GF(2^m). Jednotlive sloupce predstavuji jednotlive prvky z L a jednotlive radky jejich mocniny od 0 do t.  *)
(**)
(*     / 1        1        ... 1        \ *)
(*     | L1       L2       ... Ln       | *)
(*     | L1^2     L2^2     ... Ln^2     | *)
(* V = | .        .        .   .        | *)
(*     | .        .         .  .        | *)
(*     | .        .          . .        | *)
(*     \ L1^(t-1) L2^(t-1) ... Ln^(t-1) / *)
(**)


Unprotect[ maticeV ];
ClearAll[ maticeV ];

maticeV2[
        L_List /; ArrayDepth[L] == 2,
        M : { ired_List, p_Integer } /; ArrayDepth[ired] == 1,
        t_Integer /; t >= 1
    ] :=
    Module[ {
        ll = Length[L], matice = { },
        radek, jedna
    },
    jedna = jednotkovyPolynom[ p, Dimensions[L[[1]]] ];
    (* prvni radek matice *)
    radek = jedna& /@ Range[ll];
    AppendTo[ matice, radek ];
    Do[
        (* vynasobi radek sekvenci L *)
        radek = Thread[ krat[ radek, L, M ], List, 2 ];
        AppendTo[ matice, radek ],

        (* t-1 krat *)
        t-1
    ];
    matice
]

maticeV[
        L_List /; ArrayDepth[L] == 2,
        M : { ired_List, p_Integer } /; ArrayDepth[ ired ] == 1,
        t_Integer /; t >= 1
    ] :=
    Module[ {
        ll = Length[L], matice, radek, jedna
    },
    jedna = jednotkovyPolynom[ p, Dimensions[ L[[1]] ] ];
    (* prvni radek matice *)
    radek = jedna& /@ Range[ll];

    (* vynasobi radek sekvenci L t-1 krat *)
    NestList[
        Thread[ krat[ #, L, M ], List, 2 ]&,
        radek,
        t-1
    ]
]

Protect[ maticeV ];


(* ::Subsection::Closed:: *)
(*Matice K*)


(* ::Text:: *)
(*Matice K je soucast kontrolni matice Goppa kodu. Je to dolni trojuhelnikova matice s koeficienty ireducubilniho polynomu g.*)
(*Tato matice neni pro vypocet kontrolni matice pouzivana, ale je zde pro uplnost implementovana. *)
(**)
(*     / g_{t}   0       0       ... 0        \ *)
(*     | g_{t-1} g_{t}   0       ... 0        | *)
(*     | g_{t-2} g_{t-1} g_{t}   ... 0        | *)
(* K = | .       .       .       .   .        | *)
(*     | .       .       .        .  .        | *)
(*     | .       .       .         . .        | *)
(*     \ g_{1}   g_{2}   g_{2}   ... g_{t}    / *)
(**)


Unprotect[ maticeK ];
ClearAll[ maticeK ];

maticeK[
        g_List /; ArrayDepth[g] == 2
    ] :=
    Module[ {
        l = Length[g], nula, indexy
    },
    nula = nulovyPolynom[ 2, Dimensions[g][[2]] ];
    indexy = Table[ Max[ i-j+1, 0 ], { i, 1, l-1 }, { j, 1, l-1 } ];
    Replace[ indexy, { 0 -> nula, a_ :> g[[a]] }, {2} ]
]

Protect[ maticeK ];


(* ::Subsection::Closed:: *)
(*modifikovanyEEA*)


(* ::Text:: *)
(*Pouziva se pro dekodovani binarnich Goppa kodu. Rozdeli polynom a (modulo G) na dva polynomy Alpha a Beta tak, ze stupen(Alpha)<(stupen(G)/2) a stupen(Beta)<=stupen(G)/2*)


Unprotect[ modifikovanyEEA ];
ClearAll[ modifikovanyEEA ];

modifikovanyEEA[
        a_List /; ArrayDepth[a] == 2,
    M : { polyG_List, MVnitrni : { polyIred_List, p_Integer } } /; ArrayDepth[polyG] == 2 && ArrayDepth[polyIred] == 1
] :=
Module[ {
    (* prvky tabulky EEA *)
    A = polyG, tA = {     nulovyPolynom[ p, Rest[ Dimensions[polyG] ] ] },
    B = a,     tB = { jednotkovyPolynom[ p, Rest[ Dimensions[polyG] ] ] },
    C,         tC,
    (* pomocne promenne *)
    stupenA = stupen[polyG], stupenB = stupen[a], stupenG = stupen[polyG],
    podil
},
While[ stupenB > stupenG/2,
    { podil, C } = dlouheDeleni[ A, B, MVnitrni ];
    tC = orizniPolynom[ plus[ tA, krat[ podil, tB, M ], p] ];
    (* "posun radku" v tabulce *)
    A = B; tA = tB;
    B = C; tB = tC;
    stupenA = stupen[A];
    stupenB = stupen[B];
];
(* redukce jen pro uplnost *)
{ redukujPolynom[ B, M ], redukujPolynom[ tB, M ] }
]

Protect[ modifikovanyEEA ];


(* ::Subsection::Closed:: *)
(*invertujZakodovaniMatG*)


(* ::Text:: *)
(*Funkce vrati vektor d z opraveneho vektoru c (vektor c musi byt kodovym slovem). K "dekodovani" pouzije sloupce matice G, ktere odpovidaji dane pozici ve vektoru d ( {0,...0,1}, {0,...0,1,0}, ... ).*)


Unprotect[ invertujZakodovaniMatG ];
ClearAll[ invertujZakodovaniMatG ];

invertujZakodovaniMatG::na = "Nepodarilo se dekodovat. Vysledek (`1`) je kratsi nez k.";

invertujZakodovaniMatG[
        c_List /; ArrayDepth[c] == 1,
        G_List /; ArrayDepth[G] == 2
    ] :=
    Module[ {
        k, n, indexy, d
    },
    { k, n } = Dimensions[G];
    (* pozice samostatnych bitu v matici *)
    indexy = Flatten[ Position[ Transpose[G], #, 1 ][[1]]& /@ Reverse[Permutations[jednotkovyPolynom[ 2, k ]]]];
    d = c[[indexy]];
    If[ Length[d] != k,
        Message[ dekoduj::na, d ];
        Return[]
    ];
    d
]

Protect[ invertujZakodovaniMatG ];


(* ::Subsection::Closed:: *)
(*nahodnyChybovyVektor*)


(* ::Text:: *)
(*Nahodny vektor delky n (nad GF(2)) s Hammingovou vahou t*)


Unprotect[ nahodnyChybovyVektor ];
ClearAll[ nahodnyChybovyVektor ];

(* zamicha {1,1,..1,0,0,..,0} *)
nahodnyChybovyVektor[
        n_Integer,
        t_Integer
    ] /; t <= n :=
    RandomSample[
        Join[ 1& /@ Range[t], 0& /@ Range[n-t] ]
    ]

Protect[ nahodnyChybovyVektor ];


(* ::Subsection::Closed:: *)
(*zakoduj*)


(* ::Text:: *)
(*Zakodovani slova v libovolnem kodu -- jen maticove nasobeni (modulo 2).*)


Unprotect[ zakoduj ];
ClearAll[ zakoduj ];

zakoduj[
        a_List /; ArrayDepth[a] == 1,
        matG_List /; ArrayDepth[ matG ] == 2
    ] :=
    dotNad2[ a, matG ]

Protect[ zakoduj ];


(* ::Subsection::Closed:: *)
(*generujBinarniGoppaKod*)


(* ::Text:: *)
(*Vygeneruje a vraci parametry ireducibiliniho binarniho Goppa kodu.*)
(**)
(*Dle parametru m (vnitrni teleso GF(2^m)) a t (stupen polynomu G/pocet chyb, ktere kod opravi) vygeneruje: vnitrni konecne teleso a Goppa polynom (modul), mnozinu L (support/podpora), generujici matici a predvypocitane syndromy dle mnoziny L.*)
(**)
(*Vraci data ve formatu:*)
(*{matG,modul,podporaL,syndromyL}*)
(**)
(*Volitelny parametr verbose (True/False) urcuje, zda funkce vypise proces generovani dat. *)


Unprotect[ generujBinarniGoppaKod ];
ClearAll[ generujBinarniGoppaKod ];

generujBinarniGoppaKod::nimpl = "Neni k dispozici ireducibilni polynom g. Zkuste mensi parametry.";
generujBinarniGoppaKod::para  = "Nekompatibilni parametry: n = \!\(\*SuperscriptBox[ \(2\), \(m\) ]\) = `1`, r = t.m = `2`, k = n-r = `3`";
generujBinarniGoppaKod::polyG = "Nepodarilo se najit polynom G v `1` krocich.";
generujBinarniGoppaKod::matG  = "Nepodarilo se sestrojit matici G spravnych rozmeru ([ n, k ] = [ `1`, `2` ] != [ `3`, `4` ]).";

generujBinarniGoppaKod[
        m_Integer /; m >= 2,
        t_Integer /; t >= 2,
        verbose_Symbol : False
    ] /; BooleanQ[ verbose ] :=
    Module[ {
        p = 2, n, r, k,
        polyIred, polyG, podporaL,
        modulVnitrni, modul,
        matV, matD, (* matK, *)
        matHnadFF, matHnadGF2, matH, matG,
        polyX, syndromyL,
        i, iMax,
        tMatV, tMatD, tMatH, tMatG, tSyndromy
    },
    n = p^m; r = t*m; k = n - r;
    If[ k <= 0,
        Message[ generujBinarniGoppaKod::para, n, r, k ];
        Return[]
    ];

    modul = generujModul[ { p, m }, t ];
    If[ modul == Null,
        Message[ generujBinarniGoppaKod::nimpl ];
        Return[]
    ];
    polyG        = modul[[1]];
    modulVnitrni = modul[[2]];

    podporaL = generujNticiL[m];
    tMatV = AbsoluteTiming[
        matV = maticeV[ podporaL, modulVnitrni, t ];
    ][[1]];
    tMatD = AbsoluteTiming[
        matD = Quiet[ maticeD[ podporaL, modul ] ];
    ][[1]];
    (*matK=maticeK[polyG];*)

    tMatH = AbsoluteTiming[
        matHnadFF = dotNadFF[ matV, matD, modulVnitrni ];
        (*matHnadFF=dotNadFF[matK,matHnadFF,modulVnitrni];*)
        matHnadGF2 = Flatten[ Transpose /@ matHnadFF, 1 ];
        matH = RowReduce[ matHnadGF2, Modulus->2 ];
    ][[1]];
    tMatG = AbsoluteTiming[
        matG = NullSpace[ matH, Modulus->2 ];
    ][[1]];

    If[ verbose,
        Print[ "L = ", {podporaL } // MatrixForm ];
        (*Print["K = ",Map[MatrixForm,matK/.{nulovyPolynom[p,m]->" "},{2}]//MatrixForm];*)
        Print[ "V = ", matV // MatrixForm ];
        Print[ "D = ", Map[ MatrixForm, matD /. { nulovyPolynom[ p, m ]->" " }, {2} ] // MatrixForm ];
        (*Print["H = K.V.D = ",matHnadFF//MatrixForm];*)
        Print[ "H = V.D = ", matHnadFF // MatrixForm ];
        Print[ "H = ", matH // MatrixForm ];
        Print[ "G = ", matG // MatrixForm ];
    ];

    (* kontrola rozmeru *)
    If[ !{ k, n } === Dimensions[ matG ],
        Message[ generujBinarniGoppaKod::matG, n, k, Dimensions[matG][[2]], Dimensions[matG][[1]]];
        Return[];
    ];

    (* predpocitane dilci syndromy := {1/(x-Subscript[L, i]) mod G(x)}*)
    polyX = { jednotkovyPolynom[ p, m ], nulovyPolynom[ p, m ] };
    tSyndromy = AbsoluteTiming[
        syndromyL = inverze[ plus[ polyX, {#}, p ], modul ]& /@ podporaL;
    ][[1]];

    (* vraci nasledujici hodnoty: *)
    { { tMatV, tMatD, tMatH, tMatG, tSyndromy }, { matG, modul, podporaL, syndromyL } }
]

Protect[ generujBinarniGoppaKod ];


(* ::Subsection::Closed:: *)
(*zakodujBinarniGoppaKod*)


(* ::Text:: *)
(*Jen "alias", stejne jako pro obecny kod.*)


Unprotect[ zakodujBinarniGoppaKod ];
ClearAll[ zakodujBinarniGoppaKod ];

zakodujBinarniGoppaKod[ a_List, matG_List ] := zakoduj[ a, matG ]

Protect[ zakodujBinarniGoppaKod ];


(* ::Subsection::Closed:: *)
(*dekodujBinarniGoppaKod*)


(* ::Text:: *)
(*Pattersonuv algoritmus*)


Unprotect[ dekodujBinarniGoppaKod ];
ClearAll[ dekodujBinarniGoppaKod ];

dekodujBinarniGoppaKod[
        c_List /; ArrayDepth[c] == 1,
        GoppaKod : {
            matG_List /; ArrayDepth[ matG ] == 2,
            modul : { polyG_List, modulVnitrni : { polyIred_List, p_Integer }} /; ArrayDepth[polyG] == 2 && ArrayDepth[polyIred] == 1 && p == 2,
            podporaL_List /; ArrayDepth[ podporaL ] == 2,
            syndromyL_List /; ArrayDepth[ syndromyL ] == 3
        }
    ] /; Length[c] == Length[podporaL] == Length[syndromyL] == Length[matG[[1]]] :=
    Module[ {
        n = Length[syndromyL], m = Length[podporaL[[1]]], t = stupen[polyG],
        nula, jedna,
        polyS, polyR, polyX,
        polyAlpha, polyBeta, polySigma,
        eVektor, cOpravene, d,
        tPolyS, tPolyR, tEEA, tPolySigma, tDosazeni
    },

    nula = nulovyPolynom[ p, m ];jedna = jednotkovyPolynom[ p, m ];
    (* S(x) = Sum c_i/(x - Li) mod G(x) *)
    tPolyS = AbsoluteTiming[
        polyS = Mod[ Plus@@( c*syndromyL ), p ];
    ][[1]];
    (* pokud je syndrom nula, neni co opravovat *)
    If[ polyS === ( nula& /@ Range[t] ),

        (* then *)
        eVektor = 0& /@ Range[n],

        (* else *)
        (* R(x) = Sqrt[S(x)^-1-x] mod G(x) *)
        (* Sqrt[a] je a^((N+1)/2)kde N je pocet prvku, pro bin. telesa tedy a^2^(mt-1) *)
        polyX = { jedna, nula };
        tPolyR = AbsoluteTiming[
            polyR = umocni[ plus[ polyX, inverze[ polyS, modul ], p ], 2^(m*t-1), modul ];
        ][[1]];

        tEEA = AbsoluteTiming[
            { polyAlpha, polyBeta } = modifikovanyEEA[ polyR, modul ];
        ][[1]];

        (* Sigma(x) = Beta^2x + Alpha^2 *)
        tPolySigma = AbsoluteTiming[
            (* Alpha^2 *)
            polyAlpha = naDruhou[ #, modulVnitrni ]& /@ polyAlpha;
            polyAlpha = Riffle[ polyAlpha, {nula} ];
            (* Beta^2 *)
            polyBeta = naDruhou[ #, modulVnitrni ]& /@ polyBeta;
            polyBeta = Riffle[ polyBeta, {nula} ];
            (* Sigma(x) = Beta^2x + Alpha^2 *)
            polySigma = plus[ polyAlpha, posunPolynom[ polyBeta, 1 ], 2 ];
        ][[1]];

        (* hledani korenu *)
        (* TODO: Chien search!, nyni jen dosazeni do polynomu *)
        tDosazeni = AbsoluteTiming[
            eVektor = ( ( dosadDoPolynomu[ polySigma, #, modulVnitrni ] == nula )& /@ podporaL  ) /. { True->1, False->0 }
        ][[1]];
    ];
    cOpravene = plus[ c, eVektor, p ];

    d = invertujZakodovaniMatG[ cOpravene, matG ];
    { { tPolyS, tPolyR, tEEA, tPolySigma, tDosazeni }, { d, eVektor } }
]

Protect[ dekodujBinarniGoppaKod ];

