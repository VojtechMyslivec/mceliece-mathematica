(* ::Package:: *)

(* ::Text:: *)
(*Vojtech Myslivec, FIT CVUT v Praze, 2015/2016 *)


(* ::Title:: *)
(*Kryptosystem McEliece*)


(* ::Subsection:: *)
(*Zavislosti*)


(* ::Text:: *)
(*Nutno pridat nasledujici "baliky"*)
(*"src/rozsirenaBinarniTelesa.m"*)


(* ::Section:: *)
(*McEliece -- Implementace*)


(* ::Subsection::Closed:: *)
(*Chybove zpravy*)


McEliece::delkam = "Zprava m musi byt delky k";
McEliece::delkac = "Zprava c musi byt delky n";


(* ::Subsection::Closed:: *)
(*nahodnaPermutacniMatice*)


(* ::Text:: *)
(*Nahodna permutacni matice n*n*)


Unprotect[ nahodnaPermutacniMatice ];
ClearAll[ nahodnaPermutacniMatice ];


nahodnaPermutacniMatice[ n_ ] :=
    RandomSample[ IdentityMatrix[n] ]


Protect[ nahodnaPermutacniMatice ];


(* ::Subsection::Closed:: *)
(*nahodnaRegularniMatice*)


(* ::Text:: *)
(*Nahodna regularni matice k*k (nad GF(2)*)
(*Vraci matici M a pocet pokusu, kolik matic bylo treba vygenerovat.*)


Unprotect[ nahodnaRegularniMatice ];
ClearAll[ nahodnaRegularniMatice ];


nahodnaRegularniMatice[ k_ ] := Module[
    { r = 0, i = 0, M },
    While[ r != k,
        M = RandomChoice[ { 0, 1 }, {k, k } ];
        r = MatrixRank[ M, Modulus->2 ];
        i++
    ];
    { M, i }
]


Protect[ nahodnaRegularniMatice ];


(* ::Subsection::Closed:: *)
(*generujMcEliece*)


(* ::Text:: *)
(*Linearni kod (n,k) opravujici t chyb,*)
(* - generujici matice G*)
(* - nahodna k*k regularni matice S*)
(* - nahodna n*n permutacni matice P*)


(* ::Text:: *)
(*Parametry (n,k,t), verejny klic (hatG), soukromy klic (G, S^-1, P^-1)*)


Unprotect[ generujMcEliece ];
ClearAll[ generujMcEliece ];


generujMcEliece[
        m_Integer /; m>= 2,
        t_Integer /; t >= 2,
        verbose_Symbol : False
    ] /; BooleanQ[ verbose ] := Module[ {
        p = 2, n, k,
        GoppaKod, matG, matS, matP, hatG,
        soukromyKlic, verejnyKlic, parametry
    },
    n = p^m; k = n - m*t;

    GoppaKod = generujBinarniGoppaKod[ m, t ];
    matG = GoppaKod[[1]];
    matS = nahodnaRegularniMatice[ k ][[1]];
    matP = nahodnaPermutacniMatice[ n ];
    hatG = dotNad2[ matS, matG, matP ];

    If[ verbose == True,
        Print[
            "\!\(\*OverscriptBox[ \(G\), \(^\) ]\) = SGP = ", matS // MatrixForm, matG // MatrixForm, matP // MatrixForm ,
            "\n\!\(\*OverscriptBox[ \(G\), \(^\) ]\) =", hatG // MatrixForm
        ];
    ];

    matS = Inverse[ matS, Modulus->2 ];
    matP = Inverse[ matP, Modulus->2 ];

    soukromyKlic = { GoppaKod, matS, matP };
    verejnyKlic  = { hatG };
    parametry    = { n, k, t };

    { soukromyKlic, verejnyKlic, parametry }
]


Protect[ generujMcEliece ];


(* ::Subsection::Closed:: *)
(*sifrujMcEliece*)


(* ::Text:: *)
(*Algoritmus: "zakodovat" zpravu m (delky k), pomoci "generujici" matice*)
(* hatG a pricist nahodny chybovy vektor z (delky n) s Hammingovou vahou t*)
(* c = m hatG + z*)


Unprotect[ sifrujMcEliece ];
ClearAll[ sifrujMcEliece ];
sifrujMcEliece::delkam = McEliece::delkam;


sifrujMcEliece[
        m_List,
        verejnyKlic : { hatG_List },
        parametry : { n_Integer, k_Integer, t_Integer }
    ] := Module[ {
        p = 2, c,
        z = nahodnyChybovyVektor[ n, t]
    },
    If[ Length[ m ] != k,
        Print[ sifrujMcEliece::delkam ];
        Return[]
    ];
    c = plus[ dotNad2[ m, hatG ], z, p ];
    c
]


Protect[ sifrujMcEliece ];


(* ::Subsection::Closed:: *)
(*desifrujMcEliece*)


(* ::Text:: *)
(*Algoritmus:*)
(* - ziskat *)
(*   hatc = c P^-1 -- inverze permutace*)
(* - dekodovat *)
(*   hatm z hatc -- kod odstrani chybovy vektor*)
(* - vypocitat puvodni zpravu m=*)
(*   m = hatm S^-1*)


Unprotect[ desifrujMcEliece ];
ClearAll[ desifrujMcEliece ];
desifrujMcEliece::delkac = McEliece::delkac;


desifrujMcEliece[
        c_List,
        soukromyKlic : { GoppaKod_List, invS_List, invP_List },
        parametry : { n_Integer, k_Integer, t_Integer }
    ] := Module[ {
        hatc, hatm, m
    },
    If[ Length[ c ] != n,
        Print[ desifrujMcEliece::delkac ];
        Return[ ]
    ];
    hatc = dotNad2[ c, invP ];
    hatm = dekodujBinarniGoppaKod[ hatc, GoppaKod ][[1]];
    m = dotNad2[ hatm, invS ];
    m
]


Protect[ desifrujMcEliece ];

