(* ::Package:: *)

(* ::Text:: *)
(*Vojt\[EHacek]ch Myslivec, FIT \[CapitalCHacek]VUT v Praze, 2015/2016 *)


(* ::Title:: *)
(*Rozsirena (binarni) telesa*)


(* ::Text:: *)
(*Overovano pouze pro praci s telesy GF(q^n) s charakteristikou 2! Tedy GF(2^m^n^\[AscendingEllipsis])*)


(* ::Section:: *)
(*konecnaTelesa \[Dash] implementace*)


(* ::Subsection::Closed:: *)
(*Chybove zpravy*)


konecnaTelesa::nmale = "Delka `1` je kratsi nez (1+) stupen polynomu `2`";
konecnaTelesa::ndim = "Nekompatibilni dimenze prvku telesa a (`1`) a b (`2`).";
konecnaTelesa::nmdim = "Nekompatibilni dimenze modulu (`1`) a prvku (`2`) telesa.";
konecnaTelesa::nnula = "Polynom g nemuze byt nulovy (ani jeho nejvyssi clen).";
konecnaTelesa::nula = "Nelze spocitat inverzi k nulovemu prvku/Nelze delit nulou.";
konecnaTelesa::nimpl = "Neni implementovano nad neprvociselnymi telesy.";
konecnaTelesa::ninv = "Vysledek (`1`) neni inverze (M mozna neobsahuje ireducibilni polynom).";
konecnaTelesa::nodm = "V telese neexistuje `1`. odmocnina.";


(* ::Subsection::Closed:: *)
(*doListu*)


(* ::Text:: *)
(*Vrati seznam koeficientu polynomu p nad gf v x od nejvyssiho k nejnizsimu.*)
(*(Pripadne doplneny na n koeficientu.)*)


Unprotect[ doListu ];
ClearAll[ doListu ];
doListu::nmale = konecnaTelesa::nmale;


doListu[ a_, x_Symbol ] := Module[ {
list = Reverse[ CoefficientList[ a, x]]
},
(* osetreni prazdneho listu *)
If[ list == { },
Return[ { 0 }]
];
list
]

doListu[ a_, x_Symbol, n_Integer ] :=
Module[ { list, la },
list = CoefficientList[ a, x ];
la = Length[ list ];
If[ la>n,
Message[ doListu::nmale, n, la-1 ];
Return[ ]
];
PadLeft[ Reverse[ list ], n]
]


Protect[ doListu ];


(* ::Subsection::Closed:: *)
(*doPolynomu*)


(* ::Text:: *)
(*Vrati polynom v x dle serazeneho seznamu koeficientu v telese gf.*)


Unprotect[ doPolynomu ];
ClearAll[ doPolynomu ];
doPolynomu::nimpl = konecnaTelesa::nimpl;


doPolynomu[
koeficienty_List /; ArrayDepth[ koeficienty ] == 1,
x_Symbol
] :=
Module[ { list, n },
list = Reverse[ koeficienty ];
n = Length[ list ];
Plus@@Table[ list[[i]] x^(i-1), {i, 1, n }]
]

doPolynomu[
koeficienty_List /; ArrayDepth[ koeficienty ]>1,
x_Symbol
] :=
Message[ doPolynomu::nimpl]


Protect[ doPolynomu ];


(* ::Subsection::Closed:: *)
(*nulovyPolynom*)


(* ::Text:: *)
(*Nulovy polynom nad telesem ff "radu" r-1*)


Unprotect[ nulovyPolynom ];
ClearAll[ nulovyPolynom ];


(* ::Subsubsection::Closed:: *)
(*Pro GF(p^n)*)


(* polynom nad cisly (GF(p)). ff je konecne teleso reprezentovano cislem p *)
nulovyPolynom[ p_Integer /; p >= 2, r_Integer /; r >= 0 ] :=
0& /@ Range[ 0, r-1]


(* ::Subsubsection::Closed:: *)
(*Pro GF(q^n)*)


(* polynom nad polynomy (GF(q^n). ff je konecne teleso reprezentovano seznamem {q,n} *)
nulovyPolynom[ ff : { q_, n_Integer }, r_Integer /; r >= 0 ] :=
nulovyPolynom[ q, n ]& /@ Range[ 0, r-1]


(*
Polynom nad polynomy (GF(p^Subscript[ n, 1 ]^Subscript[ n, 2 ]^\[Ellipsis])). ff je konecne teleso reprezentovano seznamem { Subscript[ n, 1 ], Subscript[ n, 2 ], ... }
List ff odpovida hodnote Dimensions
*)
nulovyPolynom[ p_Integer /; p >= 2, ff_List ] :=
Module[ {
r = First[ ff]
},
If[ Length[ ff ] == 1,
nulovyPolynom[ p, r ],
nulovyPolynom[ p, Rest[ ff]]& /@ Range[ 1, r]
]
]


Protect[ nulovyPolynom ];


(* ::Subsection::Closed:: *)
(*jednotkovyPolynom*)


(* ::Text:: *)
(*Jednotkovy polynom nad telesem ff "radu" r-1*)


Unprotect[ jednotkovyPolynom ];
ClearAll[ jednotkovyPolynom ];


(* ::Subsubsection::Closed:: *)
(*Pro GF(p^n)*)


(* polynom nad cisly (GF(p)). ff je konecne teleso reprezentovano cislem p *)
jednotkovyPolynom[ ff_Integer /; ff >= 2, r_Integer /; r >= 0 ] :=
Append[
0& /@ Range[ 0, r-2 ],
1
]


(* ::Subsubsection::Closed:: *)
(*Pro GF(q^n)*)


(* polynom nad polynomy (GF(q^n). ff je konecne teleso reprezentovano seznamem {q,n} *)
jednotkovyPolynom[ ff : { q_, n_Integer } /; Length[ ff ] == 2, r_Integer /; r >= 0 ] :=
Append[
nulovyPolynom[ q, n ]& /@ Range[ 0, r-2 ],
jednotkovyPolynom[ q, n]
]

(* Polynom nad polynomy (GF(p^Subscript[n, 1]^Subscript[n, 2]^\[Ellipsis])). ff je konecne teleso reprezentovano seznamem {Subscript[n, 1],Subscript[n, 2],...}
List ff odpovida hodnote Dimensions *)
jednotkovyPolynom[ p_Integer /; p >= 2, ff_List ] :=
Module[ {
r = First[ ff]
},
If[ Length[ ff ] == 1,
jednotkovyPolynom[ p, r ],
Append[
nulovyPolynom[ p, Rest[ ff]]& /@ Range[ 1, r-1 ],
jednotkovyPolynom[ p, Rest[ ff]]
]
]
]


Protect[ jednotkovyPolynom ];


(* ::Subsection::Closed:: *)
(*nahodnyPolynom*)


(* ::Text:: *)
(*Nahodny polynom nad telesem ff radu (max) r-1*)


Unprotect[ nahodnyPolynom ];
ClearAll[ nahodnyPolynom ];


(* ::Subsubsection::Closed:: *)
(*Pro GF(p^n)*)


nahodnyPolynom[ ff_Integer /; ff >= 2, r_Integer /; r >= 0 ] :=
RandomInteger[ { 0, ff-1 } ]& /@ Range[ 0, r-1]


(* ::Subsubsection::Closed:: *)
(*Pro GF(q^n)*)


nahodnyPolynom[ ff : { q_, n_Integer }, r_Integer /; r >= 0 ] :=
nahodnyPolynom[ q, n ]& /@ Range[ 0, r-1]

(* Polynom nad polynomy (GF(p^Subscript[n, 1]^Subscript[n, 2]^\[Ellipsis])). ff je konecne teleso reprezentovano seznamem {Subscript[n, 1],Subscript[n, 2],...}
List ff odpovida hodnote Dimensions *)
nahodnyPolynom[ p_Integer /; p >= 2, ff_List ] :=
Module[ {
r = First[ ff]
},
If[ Length[ ff ] == 1,
nahodnyPolynom[ p, r ],
nahodnyPolynom[ p, Rest[ ff]]& /@ Range[ 1, r]
]
]


Protect[ nahodnyPolynom ];


(* ::Subsection::Closed:: *)
(*nahodnyPolynomStupne*)


(* ::Text:: *)
(*Nahodny monicky polynom nad telesem ff radu (presne) TODO r-1 (s poctem prvku r)*)


Unprotect[ nahodnyPolynomStupne ];
ClearAll[ nahodnyPolynomStupne ];


(* ::Subsubsection::Closed:: *)
(*Pro GF(p^n)*)


(* stejne jako nahodnyPolynom, jen je nutne zajistit nenulovy nejvyssi clen *)
nahodnyPolynomStupne[ ff_Integer /; ff >= 2, r_Integer /; r >= 0 ] :=
Prepend[
(* zbytek je nahodny polynom radu r-1 *)
nahodnyPolynom[ ff, r-1 ],
(* jeden jednotkovy prvek *)
1
]


(* ::Subsubsection::Closed:: *)
(*Pro GF(q^n)*)


nahodnyPolynomStupne[ ff : { q_, n_Integer }, r_Integer /; r >= 0 ] :=
Prepend[
(* zbytek je nahodny polynom radu r-1 *)
nahodnyPolynom[ ff, r-1 ],
(* jeden jednotkovy prvek *)
jednotkovyPolynom[ q, n]
]

(* Polynom nad polynomy (GF(p^Subscript[n, 1]^Subscript[n, 2]^\[Ellipsis]). ff je konecne teleso reprezentovano seznamem {Subscript[n, 1],Subscript[n, 2],...}
List ff odpovida hodnote Dimensions *)
nahodnyPolynomStupne[ p_Integer /; p >= 2, ff_List ] :=
Module[ {
r = First[ ff]
},
If[ Length[ ff ] == 1,
(* pokud je to GF(p^n), prevede na jednoduchy pripad *)
nahodnyPolynomStupne[ p, ff[[1]]],
(* jinak jeden jednotkovy prvek a zbytek nahodne *)
Prepend[
nahodnyPolynom[ p, Rest[ ff]]& /@ Range[ 1, r-1 ],
jednotkovyPolynom[ p, Rest[ ff]]
]
]
]


Protect[ nahodnyPolynomStupne ];


(* ::Subsection::Closed:: *)
(*ireducibilniPolynom*)


(* ::Text:: *)
(*Vrati ireducibilni polynom nad ff radu r-1*)


Unprotect[ ireducibilniPolynom ];
ClearAll[ ireducibilniPolynom ];
ireducibilniPolynom::nimpl = konecnaTelesa::nimpl;


(* ::Subsubsection::Closed:: *)
(*Pro GF(p^n)*)


ireducibilniPolynom[ ff_Integer /; ff >= 2, r_Integer /; r >= 2 ] :=
Module[ {
b = False, poly, x
},
While[ b == False,
poly = nahodnyPolynomStupne[ ff, r ];
b = IrreduciblePolynomialQ[ doPolynomu[ poly, x ], Modulus->ff]
];
poly
]


(* ::Subsubsection::Closed:: *)
(*Pro GF(q^n)*)


ireducibilniPolynom[ ff : { q_, n_Integer }, r_Integer /; r >= 2 ] :=
Message[ ireducibilniPolynom::nimpl]

ireducibilniPolynom[ p_Integer /; p >= 2, ff_List ] :=
If[ Length[ ff ] == 1,
ireducibilniPolynom[ p, ff[[1]]],
Message[ ireducibilniPolynom::nimpl]
]


Protect[ ireducibilniPolynom ];


(* ::Subsection::Closed:: *)
(*stupen*)


(* ::Text:: *)
(*Vrati stupen daneho polynomu. Volitelny parametr inf urcuje, zda stupen nuloveho polynomu je -\[Infinity] ci 0. *)


Unprotect[ stupen ];
ClearAll[ stupen ];


(* ::Subsubsection::Closed:: *)
(*Pro GF(q^n)*)


stupen[
a_List,
inf_ : False
] /; BooleanQ[ inf ] :=
Module[ {
nula, pozice
},
If[ ArrayDepth[ a ] == 1,
(* GF(p^n) *)
nula = 0,
(* GF(q^n), na charakteristice nezalezi *)
nula = nulovyPolynom[ 2, Rest[ Dimensions[ a]]]
];
pozice = Flatten[ Position[ a, x_ /; x != nula, {1 }]];
If[ Length[ pozice ] == 0,
(* je to nulovy polynom *)
If[ inf,
Return[ -\[Infinity]],
Return[ 0]
]
];
(* vrati stupen *)
Length[ a ]-pozice[[1]]
]


Protect[ stupen ];


(* ::Subsection::Closed:: *)
(*posunPolynom*)


(* ::Text:: *)
(*Posune polynom q o n pozic doleva. Odpovida operaci q[x]*x^n*)


Unprotect[ posunPolynom ];
ClearAll[ posunPolynom ];


posunPolynom[
p_List,
n_Integer /; n >= 0
] :=
Module[ { nula },
If[ ArrayDepth[ p ] == 1,
nula = 0,
nula = nulovyPolynom[ 2, Rest[ Dimensions[ p]]]
];
Join[ p, nula& /@ Range[ n]]
]


Protect[ posunPolynom ];


(* ::Subsection::Closed:: *)
(*orizniPolynom*)


(* ::Text:: *)
(*Orizne uvodni nuly z polynomu*)


Unprotect[ orizniPolynom ];
ClearAll[ orizniPolynom ];


orizniPolynom[
P_List
] :=
Drop[ P, Length[ P ]-stupen[ P ]-1]


Protect[ orizniPolynom ];


(* ::Subsection::Closed:: *)
(*plus*)


(* ::Text:: *)
(*Scitani polynomu v telese s charakteristikou p*)


Unprotect[ plus ];
ClearAll[ plus ];
plus::ndim = konecnaTelesa::ndim;


plus[ a_List, b_List, p_Integer /; p >= 2 ] :=
Module[ {
(* dimenze listu jsou telesa, p je charakteristika *)
la = Length[ a ], da = Dimensions[ a ],
lb = Length[ b ], db = Dimensions[ b ],
aa, bb, nula, ll
},
If[ ArrayDepth[ a ] == 1,
(* GF(p^n) *)
If[ la == lb,
Return[ Mod[ a+b, p]],
(* else -- natahni polynomy *)
ll = Max[ la, lb ];
aa = PadLeft[ a, ll ];
bb = PadLeft[ b, ll ];
Return[ Mod[ aa+bb, p]]
],

(* GF(q^n) *)
If[ !Rest[ da ] === Rest[ db ],
Message[ plus::ndim, da, db ];
Return[ ]
];
If[ la == lb,
Return[ plus[ #[[1]], #[[2]], p ]& /@ ({ a, b }\[Transpose])],
(* else -- natahni polynomy *)
nula = nulovyPolynom[ p, Rest[ da]];
ll = Max[ la, lb ];
aa = PadLeft[ a, ll, {nula } ];
bb = PadLeft[ b, ll, {nula } ];
Return[ plus[ #[[1]], #[[2]], p ]& /@ ({ aa, bb }\[Transpose])]
]
]
]


Protect[ plus ];


(* ::Subsection::Closed:: *)
(*redukujPolynom*)


(* ::Text:: *)
(*Redukce polynomu a a dle (ireducibilniho) polynomu g, nad telesem s charakteristikou p*)


Unprotect[ redukujPolynom ];
ClearAll[ redukujPolynom ];
redukujPolynom::ndim = konecnaTelesa::ndim;
redukujPolynom::nmdim = konecnaTelesa::nmdim;
redukujPolynom::nnula = konecnaTelesa::nnula;


(* ::Subsubsection::Closed:: *)
(*Pro GF(p^n)*)


(* jednoducha telesa GF(p^n) -- prevod na PolynomialMod (rychlejsi) *)
(* redukujPolynom[a,{g,p}] *)
redukujPolynom[
a_List /; ArrayDepth[ a ] == 1,
M : { g_List, p_Integer } /; ArrayDepth[ g ] == 1\[And]p >= 2
] :=
Module[ {
ll, x, ax, gx
},
(* kontrola spravneho polynomu g *)
If[
g[[1]] == 0,
Message[ redukujPolynom::nnula ];
Return[ ]
];

ll = Length[ g ]-1;(* ll je delka vysledneho prvku *)
ax = doPolynomu[ a, x ];gx = doPolynomu[ g, x ];

doListu[ PolynomialMod[ ax, {gx, p } ], x, ll]
]


(* ::Subsubsection::Closed:: *)
(*Pro GF(q^n)*)


(* M je seznam (ired.) polynomu {gVnejsi,gVnitrni} *)
redukujPolynom[
a_List /; ArrayDepth[ a ]>1,
M : { gVnejsi_List, MVnitrni_List } /; ArrayDepth[ gVnejsi ]>1\[And]Length[ MVnitrni ] == 2\[And]ArrayDepth[ gVnejsi ] == ArrayDepth[ MVnitrni[[1]]]+1
] :=
Module[ {
(* p je prvocislo/charakteristika *)
p = NestWhile[ Last, M, ListQ ], da = Dimensions[ a ], dg = Dimensions[ gVnejsi ],
gg, aa, cInv, nula, i, ll
},
(* ll je delka vysledneho prvku *)
ll = Length[ gVnejsi ]-1;
(* Kontrola spravnych rozmeru -- teles. Mohou se lisit pouze v prvni souradnici. *)
If[ !Rest[ dg ] === Rest[ da ],
Message[ redukujPolynom::nmdim, dg, da ];
Return[ ]
];

(* kontrola, ze (ired.) neni nulovy -- resp. ze jeho nejvyssi clen neni nulovy *)
nula = nulovyPolynom[ p, Rest[ dg]];
If[ gVnejsi[[1]] === nula,
Message[ redukujPolynom::nnula ];
Return[ ]
];

(* modulovany polynom -- pro jistotu *)
aa = Mod[ a, p ];
(* pokud je kratsi, neni co redukovat *)
If[ Length[ aa ] <= ll,
(* jen se (pripadne) natahne na spravnou velikost *)
Return[
PadLeft[ aa, Prepend[ Rest[ dg ], ll]]
]
];

(* prevedeni g na monicky polynom \[Dash] if kvuli zrychleni pro pripad, ze uz monicky je *)
If[ gVnejsi[[1]] === jednotkovyPolynom[ p, dg ],
(* then: je to v poradku *)
gg = gVnejsi,
(* else: je treba ho prevest na monicky *)
(* inverze nejvyssiho clenu*)
cInv = inverze[ gVnejsi[[1]], MVnitrni ];
gg = kratSkalar[ cInv, gVnejsi, MVnitrni]
];
(* natahnuty/posunuty na velikost polynomu a *)
gg = PadRight[ gg, da ];
i = 1;
While[ i <= Length[ aa ]-ll,
(* if jen mozna setri cas *)
If[ aa[[i]] != nula,
(* pricti dany nasobek (ired.) polynomu *)
aa = plus[
aa,
(* dany nasobek (skalarni) *)
(* minus pro binarni telesa neni treba *)
kratSkalar[ aa[[i]], gg, MVnitrni ],
p
]
];
gg = RotateRight[ gg ];
i++
];
Take[ aa, -ll]
]


Protect[ redukujPolynom ];


(* ::Subsection::Closed:: *)
(*krat*)


(* ::Text:: *)
(*Nasobeni polynomu v telese. Kvuli zrychleni nepouziva kratBezRedukce*)


Unprotect[ krat ];
ClearAll[ krat ];
krat::ndim = konecnaTelesa::ndim;
krat::nmdim = konecnaTelesa::nmdim;
krat::nnula = konecnaTelesa::nnula;


(* ::Subsubsection::Closed:: *)
(*Pro GF(p^n)*)


(* nasobeni v telese GF(p^n) modulo g *)
(* prevedeno na PolynomialMod -- rychlejsi *)
krat[
a_List /; ArrayDepth[ a ] == 1,
b_List /; ArrayDepth[ b ] == 1,
M : { g_List, p_Integer } /; ArrayDepth[ g ] == 1\[And]p >= 2
] :=
Module[ {
gx, ax, bx, x, ll
},
ll = Length[ g ]-1;(* delka vysledneho prvku *)
If[ g[[1]] == 0,
Message[ krat::nnula ];
Return[ ]
];

ax = doPolynomu[ a, x ];
bx = doPolynomu[ b, x ];
gx = doPolynomu[ g, x ];
doListu[ PolynomialMod[ ax bx, gx, Modulus->p ], x, ll]
]


(* ::Subsubsection::Closed:: *)
(*Pro GF(q^n)*)


(* nasobeni polynomu a a b v telese urcenem modulem M = {gVnejsi,gVnitrni} *)
krat[
a_List /; ArrayDepth[ a ]>1,
b_List /; ArrayDepth[ b ]>1,
M : { gVnejsi_List, MVnitrni_List } /; ArrayDepth[ gVnejsi ]>1\[And]Length[ MVnitrni ] == 2\[And]ArrayDepth[ gVnejsi ] == ArrayDepth[ MVnitrni[[1]]]+1
] :=
Module[ {
la = Length[ a ], da = Dimensions[ a ],
lb = Length[ b ], db = Dimensions[ b ],
p = NestWhile[ Last, M, ListQ ],
bb, i, cc, dd, nula
},
(* mohou se lisit pouze v prvni dimenzi *)
If[ !Rest[ da ] === Rest[ db ],
Message[ krat::ndim, da, db ];
Return[ ]
];
If[ !Rest[ Dimensions[ gVnejsi]] === Rest[ da ],
Message[ krat::nmdim, Dimensions[ gVnejsi ], da ];
Return[ ]
];
(* dvojnasobna dimenze -- dvojita delka pro nasobeni *)
dd = Prepend[ Rest[ da ], la+lb-1 ];
bb = PadLeft[ b, dd ];
(* v cc vznika vysledek *)
cc = nulovyPolynom[ p, dd ];
(* pomocny nulovy polynom *)
nula = nulovyPolynom[ p, Rest[ da]];
(* postupne nasobeni pres koeficienty (od nejnizsich radu \[Rule] odzadu) *)
For[ i = la, i>0, i--,
If[
!a[[i]] === nula,
(* k cc pricte nasobek bb *)
cc = plus[
cc,
kratSkalar[ a[[i]], bb, MVnitrni ],
p
]
];
bb = RotateLeft[ bb ] (* "b* = { 1, 0 }" *)
];
redukujPolynom[ cc, M]
]


Protect[ krat ];


(* ::Subsection::Closed:: *)
(*inverze*)


Unprotect[ inverze ];
ClearAll[ inverze ];
Unprotect[ inverzeITT ];
ClearAll[ inverzeITT ];
inverze::nmdim = konecnaTelesa::nmdim;
inverze::nnula = konecnaTelesa::nnula;
inverze::nula = konecnaTelesa::nula;
inverze::ninv = konecnaTelesa::ninv;


(* ::Subsubsection::Closed:: *)
(*Pro GF(p^n)*)


(* inverze v telese GF(p^n) -- prevedeno na vestaveny PolynomialExtendedGCD (je rychlejsi) *)
inverze[
a_List /; ArrayDepth[ a ] == 1,
M : { g_List, p_Integer } /; ArrayDepth[ g ] == 1\[And]p >= 2
] :=
Module[ {
(* vysledna delka prvku *)
ll = Length[ g ]-1,
gx, ax, x
},
If[ g[[1]] == 0,
Message[ inverze::nnula ];
Return[ ]
];
If[ a === nulovyPolynom[ p, Length[ a]],
Message[ inverze::nula ];
Return[ ]
];
ax = doPolynomu[ a, x ];
gx = doPolynomu[ g, x ];
doListu[ PolynomialExtendedGCD[ ax, gx, x, Modulus->p ][[2, 1]], x, ll]
]


(* ::Subsubsection::Closed:: *)
(*Pro GF(q^n)*)


(* ::Text:: *)
(*Inverze ITT*)


(* ::Text:: *)
(*Algoritmus Itoh-Teechai-Tsujii (ITT).*)
(*Vyuziva toho, ze inverze je umocnovani na 2^m-2=Subscript[(111\[Ellipsis]10), 2]  neboli: 2(2^(m-1)-1)=2Subscript[(111\[Ellipsis]1), 2] \[Dash] m-1 rika na "kolik jednicek" se ma umocnit a nakonec jeste jednou na druhou.*)
(*Vychazi z S&M algoritmu, ale ma Log(Log(n)) operaci nasobeni oproti Log(n) u S&M.*)


inverzeITT[
a_List /; ArrayDepth[ a ]>1,
M : { gVnejsi_List, MVnitrni_List } /; ArrayDepth[ gVnejsi ]>1\[And]Length[ MVnitrni ] == 2\[And]ArrayDepth[ gVnejsi ] == ArrayDepth[ MVnitrni[[1]]]+1
] :=
Module[ {
p = NestWhile[ Last, M, ListQ ],
aa = redukujPolynom[ a, M ],
inv, m, rozvoj, bb, i, k, jednotka, dim
},
If[ !Rest[ Dimensions[ gVnejsi]] === Rest[ Dimensions[ a]],
Message[ inverze::nmdim, Dimensions[ gVnejsi ], Dimensions[ a]];
Return[ ]
];
(* dimenze vysledneho prvku (jako dimenze ired. polynomu, ale o 1 mensi *)
dim = Dimensions[ aa ];
jednotka = jednotkovyPolynom[ p, dim ];
(* pokud je volano na jednotkovy polynom, vrati jednotkovy polynom *)
If[ aa === jednotka,
Return[ aa]
];
If[ aa === nulovyPolynom[ p, dim ],
Message[ inverze::nula ];
Return[ ]
];
(* rad je 2^m *)
m = Times@@dim;
(* binarni rozvoj poctu jednicek v exponentu *)
rozvoj = IntegerDigits[ m-1, 2 ];
(* v inv vznika inverze, v k je pocet hotovych jednicek v exp. *)
inv = aa;    k = 1;
(* for cyklus pres binarni rozvoj m => Log(Log(n)) iteraci *)
For[ i = 2, i <= Length[ rozvoj ], i++,
(* umocni inv (inv) na 2^k => zdvojnasobi (binarni) velikost exp. *)
bb = Nest[ naDruhou[ #, M ]&, inv, k ];
(* prinasobi puvodni hodnotu => dvojnasobna velikost exp. se samymi 1 *)
inv = krat[ bb, inv, M ];    k = 2k;
(* pokud je bit v bin. reprezentaci poctu jednicek 1 *)
If [ rozvoj[[i]] == 1,
 (* pricte jeste jednu jednicku do exponentu \[Dash] 
umocni na druhou a prinasobi puvodni a (a zvedne k) *)
inv = krat[ naDruhou[ inv, M ], aa, M ];    k = k+1;
];
];
(* posledni umocneni na druhou *)
inv = naDruhou[ inv, M ];
(* kontrola pro pripad, ze g nebyl ireducibilni polynom *)
If[ !(krat[ inv, aa, M ] === jednotka),
Message[ inverze::ninv, inv ];
Return[ ]
];
inv
]


(* ::Text:: *)
(*Inverze EEA*)


(* ::Text:: *)
(*Rozsirenym Euklidovskym algoritmem.*)


inverze[
a_List /; ArrayDepth[ a ]>1,
M : { gVnejsi_List, MVnitrni_List } /; ArrayDepth[ gVnejsi ]>1\[And]Length[ MVnitrni ] == 2\[And]ArrayDepth[ gVnejsi ] == ArrayDepth[ MVnitrni[[1]]]+1,
vypis_ : False
] /; BooleanQ[ vypis ] :=
Module[ {
p = NestWhile[ Last, M, ListQ ],
(* prvky tabulky EEA *)
A = gVnejsi,                            tA,
B = redukujPolynom[ a, M ], tB,
C,                                                 tC,
(* pomocne promenne *)
stupenA = stupen[ gVnejsi ], stupenB = stupen[ a ],
podil
},
(* vstupni kontroly *)
If[ !Rest[ Dimensions[ gVnejsi]] === Rest[ Dimensions[ a]],
Message[ inverze::nmdim, Dimensions[ gVnejsi ], Dimensions[ a]];
Return[ ]
];
If[ B === jednotkovyPolynom[ p, Dimensions[ B]],
Return[ B]
];
If[ B === nulovyPolynom[ p, Dimensions[ B]],
Message[ inverze::nula ];
Return[ ]
];

tA = {         nulovyPolynom[ p, Rest[ Dimensions[ gVnejsi]]] };
tB = { jednotkovyPolynom[ p, Rest[ Dimensions[ gVnejsi]]] };
If[ vypis, Print[ "Tabulka EEA : \n", Grid[ { {"    ", A, tA }, {"    ", B, tB }}, Frame->All, Alignment->Right, Spacings->4]]];
While[ stupenB != 0,
(* {podil,zbytek} = deleni[A/B,M] *)
{ podil, C } = dlouheDeleni[ A, B, MVnitrni ];
podil = orizniPolynom[ podil ];
C = orizniPolynom[ C ];
tC = orizniPolynom[ plus[ tA, krat[ podil, tB, M ], p]];
(* podmineny vypis *)
If[ vypis, Print[ Grid[ { {podil, C, tC }}, Frame->All, Alignment->Right, Spacings->4]]];

(* "posun radku" v tabulce *)
A = B;tA = tB;
B = C;tB = tC;
stupenA = stupen[ A ];
stupenB = stupen[ B ];
];

(* prepocet na "1" v priape, ze to bude jina jednotka okruhu *)
tC = inverze[ Last[ B ], MVnitrni ];
tB = krat[ { tC }, tB, M ];
tB
]


Protect[ inverze ];
Protect[ inverzeITT ];


(* ::Subsection::Closed:: *)
(*naDruhou*)


(* ::Text:: *)
(*Umocneni polynomu v telese na druhou.*)
(*TODO: nejspis jen pro telesa s charakteristikou 2!*)


Unprotect[ naDruhou ];
ClearAll[ naDruhou ];
naDruhou::nnula = konecnaTelesa::nnula;
naDruhou::nmdim = konecnaTelesa::nmdim;


(* ::Subsubsection::Closed:: *)
(*Pro GF(p^n)*)


(* mocneni v telese GF(p^n) modulo g *)
(* prevedeno na PolynomialMod -- rychlejsi *)
naDruhou[
a_List /; ArrayDepth[ a ] == 1,
M : { g_List, p_Integer } /; ArrayDepth[ g ] == 1\[And]p >= 2
] :=
Module[ {
gx, ax, x, ll
},
ll = Length[ g ]-1;
If[ g[[1]] == 0,
Message[ naDruhou::nnula ];
Return[ ]
];

ax = doPolynomu[ a, x ];
gx = doPolynomu[ g, x ];
doListu[ PolynomialMod[ ax ax, gx, Modulus->p ], x, ll]
]


(* ::Subsubsection::Closed:: *)
(*Pro GF(q^n)*)


(* mocneni polynomu a a b v telese urcenem modulem M = {gVnejsi,gVnitrni} *)
naDruhou[
a_List /; ArrayDepth[ a ]>1,
M : { gVnejsi_List, MVnitrni_List } /; ArrayDepth[ gVnejsi ]>1\[And]Length[ MVnitrni ] == 2\[And]ArrayDepth[ gVnejsi ] == ArrayDepth[ MVnitrni[[1]]]+1
] :=
Module[ {
p = NestWhile[ Last, M, ListQ ],
nula, cc, ll
},
If[ !Rest[ Dimensions[ gVnejsi]] === Rest[ Dimensions[ a]],
Message[ naDruhou::nmdim, Dimensions[ gVnejsi ], Dimensions[ a]];
Return[ ]
];
(* na druhou jsou puvodni koeficienty polynomu na druhou vyplnene nulami a redukovany polynomem *)
(* pomocny nulovy polynom *)
nula = nulovyPolynom[ p, Rest[ Dimensions[ gVnejsi]]];
cc = naDruhou[ #, MVnitrni ]& /@ a;
cc = Riffle[ cc, {nula } ];
redukujPolynom[ cc, M]
]


Protect[ naDruhou ];


(* ::Subsection::Closed:: *)
(*umocni*)


(* ::Text:: *)
(*Umocneni na cele cislo. Implementovano pomoci Square & Multiply.*)
(*Neni odolny vuci casovemu postranimu kanalu! Je to za cenu zrychleni*)


Unprotect[ umocni ];
ClearAll[ umocni ];
umocni::nmdim = konecnaTelesa::nmdim;



(* ::Subsubsection::Closed:: *)
(*Pro GF(q^n)*)


(* mocneni v libovolnem telese GF(q^n) *)
umocni[
a_List,
n_Integer,
M : { gVnejsi_List, MVnitrni_ } /; (ArrayDepth[ gVnejsi ] == 1\[And]ArrayDepth[ MVnitrni ] == 0)\[Or](ArrayDepth[ gVnejsi ] == ArrayDepth[ MVnitrni[[1]]]+1)
] :=
Module[ {
p = NestWhile[ Last, M, ListQ ],
aa, da, mocnina, nasobek, rozvoj, i, b1, b2
},
(* kontrola spravneho rozmeru modulu *)
If[ !Rest[ Dimensions[ gVnejsi]] === Rest[ Dimensions[ a]],
Message[ umocni::nmdim, Dimensions[ gVnejsi ], Dimensions[ a]];
Return[ ]
];

aa = redukujPolynom[ a, M ];(* redukuj \[Rule] obstara dorovnani *)
da = Dimensions[ aa ];

If[ n == 0,
Return[ jednotkovyPolynom[ p, da]]
];

If[ n>0,
(* then *)
rozvoj = IntegerDigits[ n, 2 ],
(* else *)
rozvoj = IntegerDigits[ -n, 2 ];
aa = inverze[ a, M]
];
mocnina = aa;
(* square & multiply -- zacina na druhem bitu rozvoje *)
For[ i = 2, i <= Length[ rozvoj ], i++,
(* square vzdy *)
mocnina = naDruhou[ mocnina, M ];
(* multiply, pouze pokud je odpovidajici bit 1 *)
(* MOZNY SIDE CHANNEL ATTACK! -- nutne pocitat vzdy *)
nasobek = krat[ mocnina, aa, M ];
If[ rozvoj[[i]] == 1,
mocnina = nasobek;
];
];
mocnina
]


Protect[ umocni ];


(* ::Subsection::Closed:: *)
(*odmocni*)


(* ::Text:: *)
(*Mocneni na 1/r \[Dash] pro 2. odmocninu je tak mozne pouzit odmocni[a,2,M]*)
(*TODO: nyni zkousi vsechny rady telesa, nez najde odmocninu. Prostor pro exaktni algoritmus.*)


Unprotect[ odmocni ];
ClearAll[ odmocni ];
odmocni::nodm = konecnaTelesa::nodm;


odmocni[
a_List,
r_Integer,
M : { gVnejsi_List, MVnitrni_ }
] :=
Module[ {
p = NestWhile[ Last, M, ListQ ], odmocnina,
delitele, aa, rad, invr
},
aa = redukujPolynom[ a, M ];(* redukuj \[Rule] obstara dorovnani *)

(* rad telesa (v pripade ) *)
rad = 2^Times@@Dimensions[ aa ]-1;
(* mozne rady telesa -- delitele dle Lagrange (krome 1 a radu samotneho) *)
delitele = Divisors[ rad ][[2;;-2]];

(* pokud najde odmocninu, vrati ji *)
While[ True,
(* inverze r v modulu radu *)
invr = Quiet[ PowerMod[ r, -1, rad]];
If[ IntegerQ[ invr ],
odmocnina = umocni[ aa, invr, M ];
If[ aa === umocni[ odmocnina, r, M ],
Return[ odmocnina]
];
];

(* pokud dosly moznosti, odmocnina neexistuje *)
If[ Length[ delitele ] <= 0,
Message[ odmocni::nodm, r ];
Return[ ]
];

(* jinak zkusi jiny rad *)
rad = Last[ delitele ];
delitele = Rest[ delitele ];
]
]


Protect[ odmocni ];


(* ::Subsection::Closed:: *)
(*dlouheDeleni*)


(* ::Text:: *)
(*Podil dvou polynomu P/Q s koeficienty z konecneho telesa. Pouziva algoritmus "dlouhe deleni" a vraci {podil,zbytek}.*)


Unprotect[ dlouheDeleni ];
ClearAll[ dlouheDeleni ];
dlouheDeleni::ndim = konecnaTelesa::ndim;
dlouheDeleni::nmdim = konecnaTelesa::nmdim;
dlouheDeleni::nula = konecnaTelesa::nula;


(* ::Subsubsection::Closed:: *)
(*Pro GF(p^n)*)


(* deleni polynomu nad telesem GF(p) *)
(* prevedeno na PolynomialReduce *)
dlouheDeleni[
polyP_List /; ArrayDepth[ polyP ] == 1,
polyQ_List /; ArrayDepth[ polyQ ] == 1,
p_Integer /; p >= 2
] :=
Module[ {
px, qx, x, red
},
If[ stupen[ polyQ, True ] == -\[Infinity],
Message[ dlouheDeleni::nula ];
Return[ ];
];
px = doPolynomu[ polyP, x ];
qx = doPolynomu[ polyQ, x ];
red = Flatten[ PolynomialReduce[ px, qx, x, Modulus->p]];
{ doListu[ red[[1]], x ], doListu[ red[[2]], x ] }
]


(* ::Subsubsection::Closed:: *)
(*Pro GF(q^n)*)


dlouheDeleni[
polyP_List /; ArrayDepth[ polyP ]>1,
polyQ_List /; ArrayDepth[ polyQ ]>1,
M : { gVnitrni_List, MVnitrni_ } /; (ArrayDepth[ gVnitrni ] == 1\[And]ArrayDepth[ MVnitrni ] == 0)\[Or](ArrayDepth[ gVnitrni ] == ArrayDepth[ MVnitrni[[1]]]+1)
] :=
Module[ {
p = NestWhile[ Last, M, ListQ ],
dP = Dimensions[ polyP ], dQ = Dimensions[ polyQ ],
q, koef, koefPoly,
zbytek = polyP, podil = nulovyPolynom[ 2, Dimensions[ polyQ]],
stupenZbytku = stupen[ polyP ], stupenQ = stupen[ polyQ, True ], rozdilStupnu
},
(* Vstupni kontrola *)
If[ stupenQ == -\[Infinity],
Message[ dlouheDeleni::nula ];
Return[ ]
];
If[ !Rest[ dP ] === Rest[ dQ ],
Message[ dlouheDeleni::ndim, dP, dQ ];
Return[ ]
];
If[ !Rest[ Dimensions[ gVnitrni]] === Rest[ Rest[ dP]],
Message[ dlouheDeleni::nmdim, Dimensions[ gVnitrni ], dP ];
Return[ ]
];
(* cim se bude delit *)
q = inverze[ polyQ[[-1-stupenQ]], M ];
rozdilStupnu = stupenZbytku-stupenQ;
While[ rozdilStupnu >= 0,
(* podil dvou koeficientu *)
koef = krat[ zbytek[[-1-stupenZbytku]], q, M ];
(* pricteni k celkovemu podilu *)
podil = plus[ podil, posunPolynom[ { koef }, rozdilStupnu ], p ];
(* novy zbytek *)
zbytek = plus[ zbytek, posunPolynom[ kratSkalar[ koef, polyQ, M ], rozdilStupnu ], p ];

stupenZbytku = stupen[ zbytek ];
rozdilStupnu = stupenZbytku-stupenQ;
];
(* vynechani uvodnich nul *)
zbytek = orizniPolynom[ zbytek ];
{ podil, zbytek }
]


Protect[ dlouheDeleni ];


(* ::Subsection::Closed:: *)
(*kratSkalar*)


(* ::Text:: *)
(*Nasobeni kazdeho clenu polynomu a skalarem c modulo M*)


(* ::Text:: *)
(*Neni implementovana peclive \[Dash] je to jen fce krat[] namapovana na prvky polynomu. Chyby vyhodi az vnitrni funkce krat*)


Unprotect[ kratSkalar ];
ClearAll[ kratSkalar ];


(* ::Subsubsection::Closed:: *)
(*Pro GF(p^n)*)


kratSkalar[
c_Integer,
a_List /; ArrayDepth[ a ] == 1,
p_Integer /; p >= 2
] :=
Mod[ c a, p]


(* ::Subsubsection::Closed:: *)
(*Pro GF(q^n)*)


kratSkalar[
c_List,
a_List /; ArrayDepth[ a ]>1,
M : { gVnitrni_List, MVnitrni_ } /; (ArrayDepth[ gVnitrni ] == 1\[And]ArrayDepth[ MVnitrni ] == 0)\[Or](ArrayDepth[ gVnitrni ] == ArrayDepth[ MVnitrni[[1]]]+1)
] :=
krat[ c, #, M ]& /@ a


Protect[ kratSkalar ];


(* ::Subsection::Closed:: *)
(*kratBezRedukce*)


(* ::Text:: *)
(*Nasobeni polynomu v (nad)telese*)


Unprotect[ kratBezRedukce ];
ClearAll[ kratBezRedukce ];
kratBezRedukce::ndim = konecnaTelesa::ndim;
kratBezRedukce::nmdim = konecnaTelesa::nmdim;


(* ::Subsubsection::Closed:: *)
(*Pro GF(p^n)*)


(* nasobeni v telese GF(p^n) *)
(* prevedeno na nasobeni polynomu a Mod jejich vysledku -- rychlejsi *)
kratBezRedukce[
a_List /; ArrayDepth[ a ] == 1,
b_List /; ArrayDepth[ b ] == 1,
p_Integer /; p >= 2
] :=
Module[ {
ax, bx, x
},
ax = doPolynomu[ a, x ];
bx = doPolynomu[ b, x ];
Mod[ doListu[ ax bx, x ], p]
]


(* ::Subsubsection::Closed:: *)
(*Pro GF(q^n)*)


(* nasobeni polynomu a a b v nadtelese urcenem modulem M = {gVnitrni,MVnitrnejsi} napr. {g,p} *)
kratBezRedukce[
a_List /; ArrayDepth[ a ]>1,
b_List /; ArrayDepth[ b ]>1,
M : { gVnitrni_List, MVnitrni_ } /; (ArrayDepth[ gVnitrni ] == 1\[And]ArrayDepth[ MVnitrni ] == 0)\[Or](ArrayDepth[ gVnitrni ] == ArrayDepth[ MVnitrni[[1]]]+1)
] :=
Module[ {
la = Length[ a ], da = Dimensions[ a ],
lb = Length[ b ], db = Dimensions[ b ],
p = NestWhile[ Last, M, ListQ ],
bb, i, cc, dd, nula
},
(* mohou se lisit v prvni dimenzi *)
If[ !Rest[ da ] === Rest[ db ],
Message[ kratBezRedukce::ndim, da, db ];
Return[ ]
];
If[ !Rest[ Dimensions[ gVnitrni]] === Rest[ Rest[ da]],
Message[ kratBezRedukce::nmdim, Dimensions[ gVnitrni ], da ];
Return[ ]
];

(* dvojnasobna dimenze -- dvojita delka pro nasobeni *)
dd = Prepend[ Rest[ da ], la+lb-1 ];
bb = PadLeft[ b, dd ];
(* v c vznika vysledek *)
cc = nulovyPolynom[ p, dd ];
(* pomocny nulovy polynom *)
nula = nulovyPolynom[ p, Rest[ da]];
(* postupne nasobeni pres koeficienty (od nejnizsich radu \[Rule] odzadu) *)
For[ i = la, i>0, i--,
If[
!a[[i]] === nula,
(* k c pricte nasobek bb *)
cc = plus[
cc,
kratSkalar[ a[[i]], bb, M ],
p
]
];
bb = RotateLeft[ bb ] (* "b* = { 1, 0 } alias b* = x" *)
];
cc
]


Protect[ kratBezRedukce ];


(* ::Subsection::Closed:: *)
(*dosadDoPolynomu*)


(* ::Text:: *)
(*Dosadi prvek do polynomu a vrati vysledek*)


Unprotect[ dosadDoPolynomu ];
ClearAll[ dosadDoPolynomu ];
dosadDoPolynomu::ndim = konecnaTelesa::ndim;
dosadDoPolynomu::nmdim = konecnaTelesa::nmdim;


(* ::Subsubsection::Closed:: *)
(*Pro GF(p^n)*)


dosadDoPolynomu[
polynom_List /; ArrayDepth[ polynom ] == 1,
prvek_Integer,
p_Integer /; p >= 2
] :=
Module[ { x },
Mod[ Plus[ doPolynomu[ polynom, x ] /. x->prvek ], p]
]


(* ::Subsubsection::Closed:: *)
(*Pro GF(q^n)*)


dosadDoPolynomu[
polynom_List /; ArrayDepth[ polynom ]>1,
prvek_List,
M : { gVnitrni_List, MVnitrni_ } /; (ArrayDepth[ gVnitrni ] == 1\[And]ArrayDepth[ MVnitrni ] == 0)\[Or](ArrayDepth[ gVnitrni ] == ArrayDepth[ MVnitrni[[1]]]+1)
] :=
Module[ {
dg = Dimensions[ polynom ], da = Dimensions[ prvek ], p = NestWhile[ Last, M, ListQ ],
vysledek, i
},
If[ Rest[ dg ] != da,
Message[ dosadDoPolynomu::ndim, Rest[ dg ], da ];
Return[ ]
];
If[ !Rest[ Dimensions[ gVnitrni]] === Rest[ Rest[ Dimensions[ polynom]]],
Message[ dosadDoPolynomu::nmdim, Dimensions[ gVnitrni ], Dimensions[ polynom]];
Return[ ]
];
(* Hornerovo schema *)
vysledek = polynom[[1]];
For[ i = 2, i <= Length[ polynom ], i++,
vysledek = plus[ krat[ vysledek, prvek, M ], polynom[[i]], p]
];
vysledek
]


Protect[ dosadDoPolynomu ];
