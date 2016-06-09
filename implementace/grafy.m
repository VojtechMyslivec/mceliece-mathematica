(* ::Package:: *)

(* ::Text:: *)
(*Vojtech Myslivec, FIT CVUT v Praze, 2015/2016 *)


(* ::Title:: *)
(*Mereni*)


(* ::Section:: *)
(*Vykresleni grafu*)


(* ::Text:: *)
(*Je nutne spustit notebook mereni.m ci pripadne vlozit ulozene hodnoty*)
(*(napr. ze souboru mereni/data.txt)*)
(**)
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


SetDirectory[NotebookDirectory[]];

mMin = 3; mMax = 8; tMin = 2;


(* ::Subsubsection:: *)
(*Ziskani dat pro grafy*)


m = 5;
prvek = m-mMin+1;

mereniT = mereniM[[prvek,2]];
Transpose[ { mereniT[[All,1]], mereniT[[All,2,1,1]] } ] // MatrixForm
Transpose[ { mereniT[[All,1]], mereniT[[All,2,2,1]] } ] // MatrixForm
Transpose[ { mereniT[[All,1]], mereniT[[All,2,3,1]] } ] // MatrixForm


(* ::Text:: *)
(*Tabulka hodnot: { t, generovani, sifrovani, desifrovani }*)


Transpose[ { mereniT[[All,1]], mereniT[[All,2,1,1]], mereniT[[All,2,2,1]], mereniT[[All,2,3,1]] } ] // MatrixForm


(* ::Subsubsection:: *)
(*Grafy*)


(* ::Text:: *)
(*Spolecny format pro graf*)


styl   = { Black };
osaXt  = Style[ "t",                styl ];
osaXm  = Style[ "m",                styl ];
osaYs  = Style[ "\[CHacek]as [s]",  styl ];
osaYms = Style[ "\[CHacek]as [ms]", styl ];
osaYpc = Style[ "[%]",              styl ];

moznostiZaklad = Sequence[
    ImageSize->Medium,
    AxesStyle->Directive@@styl,
    TicksStyle->Directive@@styl
];
moznostiPlot   = Sequence[
    AxesOrigin->0,
    PlotStyle->AbsolutePointSize[Medium],
    GridLines->Automatic,
    GridLinesStyle->Directive[{Dotted,Gray}],
    moznostiZaklad
];
moznostiChart  = Sequence[
    ChartLayout->"Percentile",
    ChartStyle->"DarkRainbow",
    BarSpacing->0.3,
    AxesLabel->{osaXt,osaYpc},
    moznostiZaklad
];


(* ::Text:: *)
(*Casova zavislost na parametru t pro m = 6 az 8*)


m     = Range[ 6, 8 ];
prvek = m-mMin+1;

mereniT    = mereniM[[prvek,2]];
listGen    = Reverse[ Transpose /@ Transpose[{ mereniT[[All,All,1]], mereniT[[All,All,2,1,1]] }]      ];
listSifr   = Reverse[ Transpose /@ Transpose[{ mereniT[[All,All,1]], mereniT[[All,All,2,2,1]]*1000 }] ];
listDesifr = Reverse[ Transpose /@ Transpose[{ mereniT[[All,All,1]], mereniT[[All,All,2,3,1]] }]      ];

(* format pro graf *)
legenda  = Reverse[ Style[ #, styl ]& /@ ( StringTemplate["m = `1`"][#]& /@ m )];
moznosti = Sequence[ PlotLegends->Placed[legenda,Below], moznostiPlot ];

plGen    = ListPlot[ listGen,    moznosti, AxesLabel->{ osaXt, osaYs  } ]
plSifr   = ListPlot[ listSifr,   moznosti, AxesLabel->{ osaXt, osaYms } ]
plDesifr = ListPlot[ listDesifr, moznosti, AxesLabel->{ osaXt, osaYs  } ]

soubor = StringTemplate[ "grafy/listplot_m6-8_`1`.pdf"];

Export[ soubor["generovani"],  plGen,    ImageResolution->150 ];
Export[ soubor["sifrovani"],   plSifr,   ImageResolution->150 ];
Export[ soubor["desifrovani"], plDesifr, ImageResolution->150 ];


(* ::Text:: *)
(*Casova zavislost na parametru m pro t = tMax/2*)


tPul   = Ceiling[ (Length/@mereniM[[All,2]])/2 ];
range  = Range[Length[mereniM]];
indexy = Transpose[{ range, 2 &/@ range, tPul }];
m      = range+mMin-1;

mereniT    = Part[ mereniM, #/.List->Sequence ]& /@ indexy;
listGen    = Transpose[{ m, mereniT[[All,2,1,1]]      }];
listSifr   = Transpose[{ m, mereniT[[All,2,2,1]]*1000 }];
listDesifr = Transpose[{ m, mereniT[[All,2,3,1]]      }];

barva    = RGBColor[0.880722`,0.611041`,0.142051`];
styl     = PlotStyle->Directive[{barva,AbsolutePointSize[Medium]}];
moznosti = Sequence[ PlotRange->All, styl, moznostiPlot ];

plGen    = ListPlot[ listGen,    moznosti, AxesLabel->{ osaXm, osaYs  } ]
plSifr   = ListPlot[ listSifr,   moznosti, AxesLabel->{ osaXm, osaYms } ]
plDesifr = ListPlot[ listDesifr, moznosti, AxesLabel->{ osaXm, osaYs  } ]

soubor = StringTemplate[ "grafy/listplot_tPul_`1`.pdf"];

Export[ soubor["generovani"],  plGen,    ImageResolution->150 ];
Export[ soubor["sifrovani"],   plSifr,   ImageResolution->150 ];
Export[ soubor["desifrovani"], plDesifr, ImageResolution->150 ];


(* ::Text:: *)
(*Extrapolace casove zavislosti na parametru m pro t = tMax/2*)


ClearAll[x];
f[x_]   = Fit[ listGen, {E^x}, x ];
plExtra = Plot[ f[x], {x,3,12}, AxesOrigin->{2.5,0}, AxesLabel->{osaXm,osaYs}, Evaluate[moznostiPlot] ];

plShow  = Show[ plExtra, plGen ];
Export[ "grafy/extrapolace_generovani.pdf", plShow, ImageResolution->150 ];


(* ::Text:: *)
(*Pomer casu straveny na jednotlivych operaci (pro m = 8)*)
(*Sloupcovy graf -- BarChart*)


mereniT  = mereniM[[6,2]];
hodnotyT = mereniT[[All,1]];

(* celkovy cas generovani klice *)
tGenTotal    = mereniT[[All,2,1,1]];
(* celkovy cas desifrovani *)
tDesifrTotal = mereniT[[All,2,3,1]];

(* vyznamnou cast vypoctu se stravi na *)
(* { tMatH, tSyndromy, tMatD } pri generovani (tedy na Goppa kodu) *)
vyznamneSloupce = {3,5,2};
tGenVyznamne    = mereniT[[All,2,1,2,vyznamneSloupce]];
(* { tPolyR, tDosazeni } pri desifrovani*)
vyznamneSloupce = {2,5};
tDesifrVyznamne = mereniT[[All,2,3,2,vyznamneSloupce]];

(* zbytek = celkem - soucet vyznamnych *)
tGenZbytek    = tGenTotal    - Apply[ Plus, tGenVyznamne,    {1} ];
tDesifrZbytek = tDesifrTotal - Apply[ Plus, tDesifrVyznamne, {1} ];

(* spojeni vyznamnych sloupcu a zbytku *)
tGen    = Flatten /@ Transpose[{ tGenVyznamne,    tGenZbytek    }];
tDesifr = Flatten /@ Transpose[{ tDesifrVyznamne, tDesifrZbytek }];

(* kazdy treti sloupec bude mit stitek *)
stitky  = hodnotyT/.a_:>" "/;Mod[a,5]!=0;
stitky  = ChartLabels->{stitky,None};

(* legenda pro graf *)
legendaGen    = { "matice H", "syndromy L", "matice D", "ostatn\[IAcute]" };
legendaDesifr = { "polynom R", "hled\[AAcute]n\[IAcute] ko\[RHacek]en\[URing]", "ostatn\[IAcute]" };

moznostiGen    = Sequence[ ChartLegends->Placed[legendaGen,Below],    stitky, moznostiChart ];
moznostiDesifr = Sequence[ ChartLegends->Placed[legendaDesifr,Below], stitky, moznostiChart ];

chartGen    = BarChart[ tGen,    moznostiGen    ]
chartDesifr = BarChart[ tDesifr, moznostiDesifr ]

(* nanestesti blbne vystup do pdf -- mizi popisky osy y *)
soubor = StringTemplate[ "grafy/chart_m8_`1`.png"];

Export[ soubor["generovani"],   chartGen,    ImageResolution->150 ];
Export[ soubor["desifrovani"],  chartDesifr, ImageResolution->150 ];

