(* ::Package:: *)

(* ::Text:: *)
(*Vojt\[EHacek]ch Myslivec, FIT \[CapitalCHacek]VUT v Praze, 2015/2016 *)


(* ::Title:: *)
(*Kryptosyst\[EAcute]m McEliece*)


(* ::Subsection:: *)
(*Z\[AAcute]vislosti*)


(* ::Text:: *)
(*Nutno pridat nasledujici "baliky"*)
(*"src/rozsirenaBinarniTelesa.m"*)


(* ::Section:: *)
(*McEliece \[Dash] Implementace*)


(* ::Subsection::Closed:: *)
(*Chybove zpravy*)


McEliece::delkam="Zpr\[AAcute]va m mus\[IAcute] b\[YAcute]t d\[EAcute]lky k";
McEliece::delkac="Zpr\[AAcute]va c mus\[IAcute] b\[YAcute]t d\[EAcute]lky n";


(* ::Subsection::Closed:: *)
(*nahodnaPermutacniMatice*)


(* ::Text:: *)
(*N\[AAcute]hodn\[AAcute] permuta\[CHacek]n\[IAcute] matice n*n*)


Unprotect[nahodnaPermutacniMatice];
ClearAll[nahodnaPermutacniMatice];


nahodnaPermutacniMatice[n_]:=RandomSample[IdentityMatrix[n]]


Protect[nahodnaPermutacniMatice];


(* ::Subsection::Closed:: *)
(*nahodnaRegularniMatice*)


(* ::Text:: *)
(*N\[AAcute]hodn\[AAcute] regul\[AAcute]rn\[IAcute] matice k*k (nad Subscript[\[DoubleStruckCapitalZ], 2])*)
(*Vraci matici M a pocet pokusu, kolik matic bylo treba vygenerovat.*)


Unprotect[nahodnaRegularniMatice];
ClearAll[nahodnaRegularniMatice];


nahodnaRegularniMatice[k_]:=Module[
{r=0,M,i=0},
While[
r!=k,
M=RandomChoice[{0,1},{k,k}];
r=MatrixRank[M,Modulus->2];
i++
];
{M,i}
]


Protect[nahodnaRegularniMatice];


(* ::Subsection::Closed:: *)
(*generujMcEliece*)


(* ::Text:: *)
(*Line\[AAcute]rn\[IAcute] k\[OAcute]d (n,k) opravuj\[IAcute]c\[IAcute] t chyb,*)
(* - generujici matice G*)
(* - n\[AAcute]hodn\[AAcute] k*k regul\[AAcute]rn\[IAcute] matice S*)
(* - n\[AAcute]hodn\[AAcute] n*n permuta\[CHacek]n\[IAcute] matice P*)


(* ::Text:: *)
(*Parametry (n, k, t), ve\[RHacek]ejn\[YAcute] kl\[IAcute]\[CHacek] ( *)
(*\!\(\*OverscriptBox[\(\(G\)\()\)\), \(^\)]\), soukrom\[YAcute] kl\[IAcute]\[CHacek] (G, S^-1, P^-1)*)


Unprotect[generujMcEliece];
ClearAll[generujMcEliece];


generujMcEliece[
m_Integer/;m>= 2,
t_Integer/;t>=2,
verbose_Symbol:False
]/;BooleanQ[verbose]:=Module[{
p=2,n,k,
GoppaKod,matG,matS,matP,hatG,
soukromyKlic,verejnyKlic,parametry
},
n=p^m;k=n-m t;

GoppaKod=generujBinarniGoppaKod[m,t];
matG=GoppaKod[[1]];
matS=nahodnaRegularniMatice[k][[1]];
matP=nahodnaPermutacniMatice[n];
hatG=dotNad2[matS,matG,matP];

If[verbose==True,
Print[
"\!\(\*OverscriptBox[\(G\), \(^\)]\) = SGP = ",matS//MatrixForm,matG//MatrixForm,matP//MatrixForm ,
"\n\!\(\*OverscriptBox[\(G\), \(^\)]\) =",hatG//MatrixForm
];
];

parametry={n,k,t};
verejnyKlic={hatG};
soukromyKlic={GoppaKod,Inverse[matS,Modulus->2],Inverse[matP,Modulus->2]};

{soukromyKlic,verejnyKlic,parametry}
]


Protect[generujMcEliece];


(* ::Subsection::Closed:: *)
(*sifrujMcEliece*)


(* ::Text:: *)
(*Algoritmus: "zak\[OAcute]dovat" zpr\[AAcute]vu m (d\[EAcute]lky k), pomoc\[IAcute] "generuj\[IAcute]c\[IAcute]" matice *)
(*\!\(\*OverscriptBox[\(G\), \(^\)]\) a p\[RHacek]i\[CHacek]\[IAcute]st n\[AAcute]hodn\[YAcute] chybov\[YAcute] vektor z (d\[EAcute]lky n) s Hammingovou vahou max. t*)
(*c=m *)
(*\!\(\*OverscriptBox[\(G\), \(^\)]\)+z*)


Unprotect[sifrujMcEliece];
ClearAll[sifrujMcEliece];
sifrujMcEliece::delkam=McEliece::delkam;


sifrujMcEliece[
m_List,
verejnyKlic:{hatG_List},
parametry:{n_Integer,k_Integer,t_Integer}
]:=Module[{
p=2,c,
z=nahodnyChybovyVektor[n,t]
},
If[Length[m]!=k,
Print[sifrujMcEliece::delkam];
Return[]
];
c=plus[dotNad2[m,hatG],z,p];
c
]


Protect[sifrujMcEliece];


(* ::Subsection::Closed:: *)
(*desifrujMcEliece*)


(* ::Text:: *)
(*Algoritmus:*)
(* - z\[IAcute]skat *)
(*\!\(\*OverscriptBox[\(c\), \(^\)]\)=c P^-1 \[Dash] inverze permutace*)
(* - dek\[OAcute]dovat *)
(*\!\(\*OverscriptBox[\(m\), \(^\)]\) z *)
(*\!\(\*OverscriptBox[\(c\), \(^\)]\) \[Dash] k\[OAcute]d odstran\[IAcute] chybov\[YAcute] vektor*)
(* - vypo\[CHacek]\[IAcute]tat p\[URing]vodn\[IAcute] zpr\[AAcute]vu m=*)
(*\!\(\*OverscriptBox[\(m\), \(^\)]\) S^-1*)


Unprotect[desifrujMcEliece];
ClearAll[desifrujMcEliece];
desifrujMcEliece::delkac=McEliece::delkac;


desifrujMcEliece[
c_List,
soukromyKlic:{GoppaKod_List,invS_List,invP_List},
parametry:{n_Integer,k_Integer,t_Integer}
]:=Module[{
hatc,hatm,m
},
If[Length[c]!=n,
Print[desifrujMcEliece::delkac];
Return[]
];
hatc=dotNad2[c,invP];
hatm=dekodujBinarniGoppaKod[hatc,GoppaKod][[1]];
m=dotNad2[hatm,invS];
m
]


Protect[desifrujMcEliece];
