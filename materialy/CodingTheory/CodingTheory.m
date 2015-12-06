(*:Title:Coding Theory*)

(* :Author: Gachkov Igor, Karlstad  University Sweden *)
(* :Mathematica Version: 4.2 *) 
(* :Package Version: 1.0 *) 
(* :Context: CodingTheory`
*)
(* :Summary:
   The package consists of two parts: 
   one part with illustrative explanations, 
   the illustrative part  (commands starting with Show...)
   is considered to visualize the theoretical aspects of encoding / decoding,
   construct shift-register circuits etc,
   and another one for scientific purposes.
   Information and description:
   http://www.ingvet.kau.se/~igor/
*)          
(* :History:
   The package in MATHEMATICA in the field "Coding Theory"
    was developed for course " Error-Correcting codes with MATHEMATICA "
*)
(* Any comments or bug reports should be forwarded to one of the following:

   Igor Gashkov
Karlstad University, 
Department of Engineering Sciences, 
Physics and Mathematics 
65188  Karlstad Sweden

Igor.Gachkov@kau.se
 
+46(0)54-700-11-03
*)
		(**Coding Theory**)

BeginPackage["CodingTheory`"]

Print[Date[][[1]],"-",Date[][[2]],"-",Date[][[3]]]
Print["Igor Gachkov Karlstad University Sweden"] 
Print["Now loading ..."]

Unprotect[
AugmentedCode, 
BinaryGaloisField, 
BinaryIrreduciblePolynomials, 
BurstRectangularCodeVector, 
CyclicCode, 
DecHammingCode, 
DecRectangularCode, 
DimensionCode, 
DistanceCode, 
DualCode, 
EquivalentCode, 
ErrorTrappingDecodeBCHCode, 
ErrorTrappingDecodeHammingCode, 
ExpurgatingCode, 
ExtendedCode, 
GeneratorMatrix, 
GeneratorPolynomials, 
GilbertVarshamovBound, 
HammingBound, 
HammingCode, 
HammingMatrix, 
McWilliamsIdentity, 
MeggittDecode, 
NonsystematicEncodeCyclicCode, 
ParityCheckMatrix, 
PuncturedCode, 
RectangularCode, 
ShowBinaryGaloisField, 
ShowBurstRectangularCodeVector, 
ShowCorrectBurstRectangularCode, 
ShowCyclicCode, 
ShowDecHammingCode, 
ShowDecRectangularCode, 
ShowErrorTrappingDecoderBCHCode, 
ShowHammingCode, 
ShowMeggittDecoder, 
ShowNonsystematicEncoderCyclicCode, 
ShowReconstructInformRectangularCode, 
ShowRectangularCode, 
ShowStandardArrays, 
ShowSystematicEncode, 
ShowSystematicEncoderCyclicCode, 
ShowTableMult, 
ShowTablePlus, 
SingletonBound, 
Syndrome, 
SystematicEncode, 
SystematicEncodeCyclicCode, 
VEC, 
WeightPolynomial
]


(*AugmentedCode*)

AugmentedCode::usage =
	"AugmentedCode[H] gives the parity check matrix 
of the code augmented by adding new codewords 
( by adding the all-ones vector )."

(*BinaryGaloisField*)

BinaryGaloisField::usage =
	"BinaryGaloisField[g,x,b,a] gives the lists of the elements of
the Galois Field which has been constructed using 
the irreducible polynomial g(x), 
with primitive element b , minimal polynomials."

(*BinaryIrreduciblePolynomials*)

BinaryIrreduciblePolynomials::usage =
	"BinaryIrreduciblePolynomials[n,x] gives the list of the binary
irreducible polynomials  r(x) of degree n over the finite field 
GF(2) ={0,1}."

(*BurstRectangularCodeVector*)

BurstRectangularCodeVector::usage =
	"BurstRectangularCodeVector[n,m,inf] encodes
of the information word inf into 
the n*m rectangular code using diagonalewize 
represantation of the code words."

(*CyclicCode*)

CyclicCode::usage =
	"CyclicCode[g,n,x] gives the list { H, G, h(x) } , where
H = the parity check matrix, G = the generator matrix, 
h(x) = the parity check polynomial of the cyclic code 
with generator polynomial g and length n."

(*DecHammingCode*)

DecHammingCode::usage =
	"DecHammingCode[m,w] decodes 
the received word w and reconstructs the information word, gives the list {v,inf}"

(*DecRectangularCode*)

DecRectangularCode::usage =
	"DecRectangularCode[n,m,w] decodes the received word w
(if it is possible) and reconstructs the information word of 
the n*m rectangular code ."

(*DimensionCode*)

DimensionCode::usage =
	"DimensionCode[H] gives the dimension of the code with
parity check matrix H."

(*DistanceCode*)

DistanceCode::usage =
	"DistanceCode[H] gives the code distance of the code with
parity check matrix H."

(*DualCode*)

DualCode::usage =
	"DualCode[H] gives the parity check matrix of the code."

(*EquivalentCode*)

EquivalentCode::usage =
	"EquivalentCode[H,q] gives the parity check matrix of the equivalent 
code, constructed using the permutation q={k1, k2, ..., kn} where ki = 1-n ."

(*ErrorTrappingDecodeBCHCode*)

ErrorTrappingDecodeBCHCode::usage =
	"ErrorTrappingDecodeBCHCode[w] gives the method and 
shows the trace of the received word w into the error trapping  decoder 
of the BCH code of length 15 with generator polynomial 1+x^4+x^6+x^7+x^8. "

(*ErrorTrappingDecodeHammingCode*)

ErrorTrappingDecodeHammingCode::usage =
	"ErrorTrappingDecodeHammingCode[w,x] gives the method and 
show the trace of the received word w into the error trapping  decoder 
of the Hamming code of length 2^m  where m = Log[2,Length[w]+1] with generator polynomial 
(m == 2) g = 1 + x + x^2, 
(m == 3) g = 1 + x + x^3, 
(m == 4) g = 1 + x + x^4, 
(m == 5) g = 1 + x^2 + x^5, 
(m == 6) g = 1 + x + x^6, 
(m == 7) g = 1 + x + x^7."

(*ExpurgatingCode*)

ExpurgatingCode::usage =
	"ExpurgatingCode[H] gives the parity check matrix of a new code. 
(We expurgate K by throwing away the codewords of odd weight )."

(*ExtendedCode*)

ExtendedCode::usage =
	"ExtendedCode[H] gives the parity check matrix of the extended code."

(*GeneratorMatrix*)

GeneratorMatrix::usage =
	"GeneratorMatrix[H] gives the generator matrix of the code with
parity check matrix H."

(*GeneratorPolynomials*)

GeneratorPolynomials::usage =
	"GeneratorPolynomials[n,x] gives the list of the  generator 
polynomials of the cyclic code of length n."

(*GilbertVarshamovBound*)

GilbertVarshamovBound::usage =
	"GilbertVarshamovBound[n,k,d] gives the Gilbert - Varshamov bound  
( lower bound on the size of the code ) . 
n - length of the code, k - dimension of the code , 
d - code distance of the code 
( n,k,d - numbers or symbols )."

(*HammingBound*)

HammingBound::usage =
	"HammingBound[n,k,d] gives the Hamming bound  
( upper bound on the size of the code ) . 
n - length of the code, k - dimension of the code , 
d - code distance of the code 
( n,k,d - numbers or symbols )."

(*HammingCode*)

HammingCode::usage =
	"HammingCode[m,inf] encodes the
information inf into a Hamming code if length 2^m - 1."

(*HammingMatrix*)

HammingMatrix::usage =
	"HammingMatrix[m] gives the parity check matrix of 
the Hamming code of length 2^m - 1."

(*McWilliamsIdentity*)

McWilliamsIdentity::usage =
	"McWilliamsIdentity[f,n,x] gives the weight polynomial of 
the dual code for the code K with weight polynomial f 
and  length n, using the Mac Williams Identity ."

(*MeggittDecode*)

MeggittDecode::usage =
	"MeggittDecode[g,w,n,x] gives the method and 
show the trace of the received word w into the Meggitt decoder 
of the cyclic code with generator polynomial g(x) and length n."

(*NonsystematicEncodeCyclicCode*)

NonsystematicEncodeCyclicCode::usage =
	"NonsystematicEncodeCyclicCode[g,inf,n,x] gives the method and shows
the trace of the information word inf into the nonsystematic encoder 
of the cyclic code with generator polynomial g(x) 
and length n."

(*ParityCheckMatrix*)

ParityCheckMatrix::usage =
	"ParityCheckMatrix[G] gives the parity check matrix of the code with
generator matrix G."

(*PuncturedCode*)

PuncturedCode::usage =
	"PuncturedCode[H,n] gives the parity check matrix 
of the puncturing of a code by deleting the n-th  coordinate."

(*RectangularCode*)

RectangularCode::usage =
	"RectangularCode[n,m,inf]  encodes information into
a n*m rectangular code ."

(*ShowBinaryGaloisField*)

ShowBinaryGaloisField::usage =
	"ShowBinaryGaloisField[g,x,b,a] shows the lists of the elements of
the Galois Field which has been constructed using 
the irreducible polynomial g(x), 
with primitive element b , minimal polynomials."

(*ShowBurstRectangularCodeVector*)

ShowBurstRectangularCodeVector::usage =
	"ShowBurstRectangularCodeVector[n,m,inf,X] shows the method  
of encoding of the information word inf into 
the n*m rectangular code using diagonalewize 
represantation of the code words."

(*ShowCorrectBurstRectangularCode*)

ShowCorrectBurstRectangularCode::usage =
	"ShowCorrectBurstRectangularCode[n,m,v,X] shows the method  
of correcting burst errors of the received word w using 
the n*m rectangular code and diagonalewize 
represantation of the code words."

(*ShowCyclicCode*)

ShowCyclicCode::usage =
	"ShowCyclicCode[g,n,x] shows the method of construction of the parity
check matrix, generator matrix and parity check polynomial of the cyclic code
with generator polynomial g and length n."

(*ShowDecHammingCode*)

ShowDecHammingCode::usage =
	"ShowDecHammingCode[m,w] shows the method of decoding 
of the received word w and reconstructs the information word."

(*ShowDecRectangularCode*)

ShowDecRectangularCode::usage =
	"ShowDecRectangularCode[n,m,w] shows the method of decoding of
the received word w (if it is possible) and reconstructs the information 
word of the n*m rectangular code ."

(*ShowErrorTrappingDecoderBCHCode*)

ShowErrorTrappingDecoderBCHCode::usage =
	"ShowErrorTrappingDecoderBCHCode constructs the error trapping 
decoder of the BCH code of length 15 with 
the generator polynomial 1+x^4+x^6+x^7+x^8."

(*ShowHammingCode*)

ShowHammingCode::usage =
	"ShowHammingCode[m,inf] shows the method of encoding 
information into a Hamming code of length 2^m - 1."

(*ShowMeggittDecoder*)

ShowMeggittDecoder::usage =
	"ShowMeggittDecoder[g,n,x] constructs the 
Meggitt decoder of the cyclic code with generator polynomial g(x) 
and length n."

(*ShowNonsystematicEncoderCyclicCode*)

ShowNonsystematicEncoderCyclicCode::usage =
	"ShowNonsystematicEncoderCyclicCode[g,n,x] constructs the 
nonsystematic encoder of the cyclic code with generator polynomial g(x) 
and length n."

(*ShowReconstructInformRectangularCode*)

ShowReconstructInformRectangularCode::usage =
	"ShowReconstructInformRectangular[n,m,w,X] shows the method  
of reconstructing the information word ( where v is the corrected code word )
using the n*m rectangular code and diagonalewize 
represantation of the code words."

(*ShowRectangularCode*)

ShowRectangularCode::usage =
	"ShowRectangularCode[n,m,inf,X] shows the method of encoding information
into a n*m rectangular code ."

(*ShowStandardArrays*)

ShowStandardArrays::usage =
	"ShowStandardArrays[H,w] shows the method of decoding of
the received word w  using the Standard Array 
( H - parity check matrix of code ) ."

(*ShowSystematicEncode*)

ShowSystematicEncode::usage =
	"ShowSystematicEncode[H,inf,x] shows the method of encoding information
( if it is possible) into a systematic codeword (H - parity check matrix) ."

(*ShowSystematicEncoderCyclicCode*)

ShowSystematicEncoderCyclicCode::usage =
	"ShowSystematicEncoderCyclicCode[g,n,x] constructs the 
systematic encoder of the cyclic code with generator polynomial g(x) 
and length n."

(*ShowTableMult*)

ShowTableMult::usage =
	"ShowTableMult[n] shows the table of multiplication by modulo n."

(*ShowTablePlus*)

ShowTablePlus::usage =
	"ShowTablePlus[n] shows the table of addition by modulo n."

(*SingletonBound*)

SingletonBound::usage =
	"SingletonBound[n,k,d] gives the Singleton bound 
( upper bound on the size of the code ) . 
n - length of the code, k - dimension of the code ,
d - code distance of the code 
( n,k,d - numbers or symbols )."

(*Syndrome*)

Syndrome::usage =
	"Syndrome[H,v] gives the syndrome of the recieved word v for the code
with parity check matrix H."

(*SystematicEncode*)

SystematicEncode::usage =
	"SystematicEncode[H,inf,x] encodes information
( if it is possible) into a systematic codeword (H - parity check matrix) ."

(*SystematicEncodeCyclicCode*)

SystematicEncodeCyclicCode::usage =
	"SystematicEncodeCyclicCode[g,inf,n,x] gives the method and 
show the trace of the information word inf into the systematic encoder 
of the cyclic code with generator polynomial g(x) and length n."

(*VEC*)

VEC::usage =
	"VEC[n,m] gives the binary vector of length m ,
as the binary representation of the number n."

(*WeightPolynomial*)

WeightPolynomial::usage =
	"WeightPolynomial[H,x] gives the weight polynomial of the code with
parity check matrix H ."

Begin["`Private`"]

AugmentedCode[M_List]:=
ParityCheckMatrix[Prepend[GeneratorMatrix[M],
Table[1,{Length[M[[1]]]}]]]

BinaryGaloisField[g_,x_Symbol,b_Symbol,a_Symbol] := 
  Module[{n, o, od, br, pr, W, KL, j, i, prv, r, 
    f, v, len, minpol, probel, FPOL,Pole},Pole={};
   n = Exponent[g, x]; 
    If[n == 1, Print["Cive the another\
        polynomial with degree more then 1"]; 
      Break[]]; FL = FactorList[g, Modulus -> 2]; 
    If[Length[FL] != 2 || FL[[2,2]] != 1, 
     Print["Polynomial ", g, 
       " is not irreducible polynomial"]; Break[]
      ]; o = {}; 
Do[AppendTo[o,0], {n}]; 
    od = Take[Prepend[o,1], n]; 
    br = 0; pr = " "; W = 0; 
    Do[If[br == 0, v = Reverse[VEC[i, n]]; 
       f = Sum[v[[j]]*x^(j - 1), {j, n}]; 
       Do[If[Mod[2^n - 1, j] == 0, 
         If[j != 2^n - 1, 
           r = 
            PolynomialMod[f^j, g, 
             Modulus -> 2]]; 
          If[r == 1 && j < 2^n - 1, Break[], 
           pr = f; br = 1]], {j, 2^n - 1}]], 
     {i, 2^n - 1}]; Do[f = 
       PolynomialMod[pr^i, g, Modulus -> 2] /. 
        x -> a; FPOL = 
       FactorList[x^(2^n - 1) - 1, Modulus -> 2]\
       ; stop = 0; Do[If[stop == 1, Break[]]; 
        ar = 
         AlgebraicRules[(g /. x -> a) == 0, 
          {x, a}]; 
        If[PolynomialMod[PolynomialMod[FP\
               OL[[j,1]] /. x -> f, 2] /. ar\
            , 2] == 0, 
         minpol = FPOL[[j,1]]; stop = 1; 
          Break[]], {j, Length[FPOL]}]; 
      KL = CoefficientList[f, {a}]; prv = {}; 
      Do[If[KL[[i]] == 0, 
        AppendTo[prv,0], 
       AppendTo[prv,1]], 
       {i, Length[KL]}]; 
      len =Length[prv]; 
      If[len < n, Do[AppendTo[prv,0], {n - len}]]\
       ;If[W==0,AppendTo[Pole,{"-oo", o, 0, 0, x}];
     AppendTo[Pole,{0,od,1,1,1 + x}]];
 W = 1;AppendTo[Pole,{i,prv,b^i,f,minpol}], 
     {i, 1, 2^n - 2}];Pole]

BinaryIrreduciblePolynomials[n_Integer,x_Symbol] := 
  Module[{k, v, f, j, i, FL,BIP}, BIP = {}; 
    Do[v = VEC[i, n + 1]; If[Mod[Apply[Plus, v], 2] != 0&&Last[v]==1, 
       f = Sum[v[[j]]*x^(j - 1), {j, n + 1}]; 
         FL = FactorList[f, Modulus -> 2]; 
          If[Length[FL] == 2 && FL[[2,2]] == 1,AppendTo[BIP, f]]], 
     {i, 1, 2^(n + 1) - 1, 2}];BIP]

BurstRectangularCodeVector[n_Integer, m_Integer, inf_List] := 
  Module[{a, i, j, w, b, c}, If[TEST2[inf] == 1, 
     Print["CHECK VECTOR inf ON 0 OR 1"]; Break[]]; a = inf; 
    If[Length[a] != n*m, Print["THE LENGTH INF =!= n*m"]; Break[]]; 
   Do[a = Insert[a, X, i*(n + 1)], {i, m}]; 
    Do[AppendTo[a, X], {i, n + 1}]; a = Partition[a, n + 1];   
    Do[a[[i,n + 1]] = sum2[Drop[a[[i]], {n + 1}]], {i, m}]; 
    Do[a[[m + 1,j]] = sum2[Table[a[[i,j]], {i, m}]], {j, n + 1}]; 
     w = Table[0, {(n + 1)*(m + 1)}]; 
    Do[Do[w[[Mod[(m + 1)*(j - i) + i - 1, (m + 1)*(n + 1)] + 1]] = a[[i,j]], 
      {i, m + 1}], {j, n + 1}]; w]

CyclicCode[g_, n_Integer,x_Symbol] := 
  Module[{g1, y, Y, l, k, K, G, i, j, R, H,HM,GM,h},   
    g1 = PolynomialMod[g, 2];
    l = Length[CoefficientList[PolynomialMod[x^n - 1, g1, 
        Modulus -> 2], {x}]]; 
    If[l != 0, Print["The polynomial g =", g1, 
       " isn`t a generator polynomial, as"]; 
      Print["the generator polynomial of cyclic code of
length ", n, " must divides  ", x^n - 1]; Break[]];    
    k = n - Length[CoefficientList[g1, {x}]] + 1; 
    GM = Array[a, {k, n}]; 
    Do[L = CoefficientList[g1*x^(i - 1), {x}]; 
      Do[AppendTo[L, 0], {k - i}]; GM[[i]] = L, {i, k}];      
    h = PolynomialMod[PolynomialQuotient
[x^n - 1, g1, x],2];
    HM = Array[b, {n - k, n}]; 
    Do[R = Reverse[CoefficientList[h, {x}]]; 
      Do[AppendTo[R, 0], {n - k - i}]; 
      Do[PrependTo[R, 0], {i - 1}]; HM[[n - k - i + 1]] = R
      , {i, n - k}];{HM,GM,h}]

DecHammingCode[m_Integer, w_List] := 
  Module[{M, NS, ER, S}, 
   If[TEST2[w] == 1, Print["Check w on 0 or 1"]; Break[]]; 
    If[Length[w] != 2^m - 1, 
     Print["Length w = ", Length[w], " =!= ", 2^m - 1]; 
      Break[]];  M = HammingMatrix[m]; 
    S = Mod[M . w, 2]; 
    NS = NullSpace[M, Modulus -> 2]; 
    ER = Sum[S[[i]]*2^(m - i), {i, m}]; 
    If[ER =!= 0,W = w; W[[ER]] = Mod[W[[ER]] + 1, 2];       
     {W,LinearSolve[Transpose[NS], W, 
        Modulus -> 2]}, W = w;      
    {W,LinearSolve[Transpose[NS], W, 
        Modulus -> 2]}]]

DecRectangularCode[n_Integer, m_Integer, w_List] := 
  Module[{s, A, d, OX, OY, er, i, j, DW, a, b, p}, 
   T = 0; If[TEST2[w] == 1, 
     Print["Check w on 0 or 1"]; Break[]]; 
    If[Length[w] != (n + 1)*(m + 1), 
     Print["Check the length of w"]; Break[]];     
    s = Partition[w, n + 1];    
    OX = {}; OY = {}; d = 1; 
    Do[er = sum2[s[[i]]]; 
      If[er =!= 0, 
       AppendTo[OX, i]], {i, m + 1}]; 
    Do[er = sum2[Table[s[[i,j]], {i, m + 1}]]; 
      If[er =!= 0,  
       AppendTo[OY, j]], {j, n + 1}];
    If[Length[OX] > 1 || Length[OY] > 1, 
     Print["2 or more errors."], 
     If[Length[OX] == 0, DW = s; d = 0, 
      p = s[[OX[[1]],OY[[1]]]]; 
       s[[OX[[1]],OY[[1]]]] = Mod[p + 1, 2]; DW = s; d = 0]
      ];If[d == 1, Break[], 
     a = Drop[DW, {m + 1}]; b = Transpose[a]; 
      a = Drop[b, {n + 1}]; b = Transpose[a]; 
      a = Flatten[b]]]

DimensionCode[H_List]:=Length[GeneratorMatrix[H]]

DistanceCode[M_List] := 
  Module[{a, i, d, m, n, k, b},
If[TEST1[M] == 1, Print["Check the matrix M"]; Break[]]; 
    If[TEST2[M] == 1, Print["Check the matrix on 0 or 1"]; 
      Break[]];m=DimensionCode[M];n=Length[M[[1]]];
If[m>n-m,
Exponent[
McWilliamsIdentity[
WeightPolynomial[
DualCode[M],x],n,x][[2]],x],d = n; 
    Do[a = VEC[i, m]; k = 
       Apply[Plus,Mod[a.GeneratorMatrix[M], 2]];
If[k<d&&0<k,d=k], {i, 0, 2^m - 1}];d]]

DualCode[M_List]:=NullSpace[M,Modulus->2]

EquivalentCode[M_List,r_List]:=
Transpose[Table[Transpose[M][[r[[i]]]],{i,1,Length[r]}]]

ErrorTrappingDecodeBCHCode[W_List] := 
  Module[{P, gr1, gr2, gr3, gr4, gr5, gr6, gr7, gr8, gr9, gr10, gr11, gr122, 
      gr133, gr144, gr155, grtext, j, i, rw, correct, cor, crt, WCOR, wc, DWD,
       dw, Pnew, pr, con, cr}, w = Reverse[W]; 
    If[Length[w] != 15, 
      Print["The length of w = ", Length[w], " and != 15."];
      Break[]]; 
    If[TEST2[w] == 1, Print["Check the word w on 0 or 1."]; Break[]];
    rw = ""; correct = 0; cor = 0; crt = 0;
    WCOR = w;
    wc = ""; DWD = w; dw = "";
    Do[If[w[[i]] == 0, rw = StringInsert[rw, "0", -1], 
        rw = StringInsert[rw, "1", -1]], {i, 15}];
    P = Table["", {8}]; Pnew = P;
    gr3 = Table[Graphics[{}], {3}];
    graf = Table[Graphics[{}], {i, 3}, {j, 15}];
    gr7 = Table[Graphics[{}], {2}]; gr8 = gr7;
    Do[Do[
        gr1 = Graphics[{AbsoluteThickness[1 + If[i == 1, w[[16 - j]], 0]], 
              Line[{{-2, 7}, {63, 7}, {63, 3}}], 
              Line[{{0, 
                    7}, {0, -19}, {19, -19}, {19, -22}, {59, -22}, {59, -18},
{19, -18}, {19, -20}}]}];
        Pnew[[1]] = 
          Mod[If[i == 1, w[[16 - j]], 0] + If[IntegerQ[P[[8]]], P[[8]], 0] + 
              correct, 2];
        gr2 = 
          Graphics[{AbsoluteThickness[1 + Pnew[[1]]], Circle[{63, 2}, 1], 
              Text["+", {63, 2}], 
              Line[{{62, 2}, {1, 
                    2}, {1, -2}, {1, -2}, {3, -2}, {3, -1}, {5, -1}, {5, -3}, 
{3, -3}, {3, -1}}], Line[{{32, 2}, {32, -1}}], Line[{{48, 2}, {48, -1}}], 
              Line[{{56, 2}, {56, -1}}], 
              If[i == 1, Line[{{4, -3}, {4, -5}, {3, -6}}], 
                Line[{{4, -3}, {4, -8}}]]}];
        Do[Pnew[[l + 1]] = P[[l]];
          
          gr3[[l]] = 
            Graphics[{AbsoluteThickness[
                  1 + If[IntegerQ[Pnew[[l + 1]]], Pnew[[l + 1]], 0]], 
                Line[{{5 + 8*(l - 1), -2}, {11 + 8*(l - 1), -2}, {11 + 
                        8*(l - 1), -1}, {13 + 8*(l - 1), -1}, {13 + 
                        8*(l - 1), -3}, {11 + 8*(l - 1), -3}, {11 + 
                        8*(l - 1), -1}}], 
                If[i == 1, 
                  Line[{{12 + 8*(l - 1), -3}, {12 + 8*(l - 1), -5}, {11 + 
                          8*(l - 1), -6}}], 
                  Line[{{12 + 8*(l - 1), -3}, {12 + 8*(l - 1), -8}}]]}], {l, 
            3}];
        gr4 = 
          Graphics[{AbsoluteThickness[1 + If[IntegerQ[P[[4]]], P[[4]], 0]], 
              Line[{{29, -2}, {31, -2}}]}];
        Pnew[[5]] = Mod[Pnew[[1]] + If[IntegerQ[P[[4]]], P[[4]], 0], 2]; 
        gr5 = Graphics[{AbsoluteThickness[
                1 + Mod[Pnew[[1]] + If[IntegerQ[P[[4]]], P[[4]], 0], 2]], 
              Circle[{32, -2}, 1], Text["+", {32, -2}], 
              Line[{{33, -2}, {35, -2}, {35, -1}, {37, -1}, {37, -3}, {35, 
-3}, {35, -1}}], 
              If[i == 1, Line[{{36, -3}, {36, -5}, {35, -6}}], 
                Line[{{36, -3}, {36, -8}}]]}];
        Pnew[[6]] = P[[5]];
        gr6 = 
          Graphics[{AbsoluteThickness[
                1 + If[IntegerQ[Pnew[[6]]], Pnew[[6]], 0]], 
              Line[{{37, -2}, {43, -2}, {43, -1}, {45, -1}, {45, -3}, {43, 
-3}, {43, -1}}], 
              If[i == 1, Line[{{44, -3}, {44, -5}, {43, -6}}], 
                Line[{{44, -3}, {44, -8}}]]}];
        Do[
          gr7[[i]] = 
            Graphics[{AbsoluteThickness[
                  1 + If[IntegerQ[P[[i + 5]]], P[[i + 5]], 0]], 
                Line[{{45 + 8*(i - 1), -2}, {47 + 8*(i - 1), -2}}]}], {i, 
            2}];
        Do[
          Pnew[[l + 6]] = 
            Mod[Pnew[[1]] + If[IntegerQ[P[[l + 5]]], P[[l + 5]], 0], 2]; 
          gr8[[l]] = 
            Graphics[{AbsoluteThickness[1 + Pnew[[l + 6]]], 
                Circle[{48 + 8*(l - 1), -2}, 1], 
                Text["+", {48 + 8*(l - 1), -2}], 
                Line[{{49 + 8*(l - 1), -2}, {51 + 8*(l - 1), -2}, {51 + 
                        8*(l - 1), -1}, {53 + 8*(l - 1), -1}, {53 + 
                        8*(l - 1), -3}, {51 + 8*(l - 1), -3}, {51 + 
                        8*(l - 1), -1}}], 
                If[l == 1, 
                  If[i == 1, Line[{{52, -3}, {52, -5}, {51, -6}}], 
                    Line[{{52, -3}, {52, -8}}]], 
                  Line[{{60, -3}, {60, -15}, {62, -15}}]]}], {l, 2}]; 
        gr9 = Graphics[{AbsoluteThickness[
                1 + If[IntegerQ[P[[8]]], P[[8]], 0]], 
              Line[{{61, -2}, {63, -2}, {63, 1}}]}];
        cor = 
          If[i > 1, 
            If[Take[Pnew, 7] == {0, 0, 0, 0, 0, 0, 0} || 
                Apply[Plus, Take[Pnew, 7]] == 1, 1, 0], 0];
        gr10 = 
          Graphics[{AbsoluteThickness[1 + cor], 
              Line[{{1, -8}, {55, -8}, {55, -14}, {1, -14}, {1, -8}}], 
              Line[{{28, -14}, {28, -17}, {62, -17}}]}];
        correct = If[cor == 1 && Pnew[[8]] == 1, 1, 0];
        gr11 = 
          Graphics[{RGBColor[crt, 0, 0], AbsoluteThickness[1 + crt], 
              Line[{{62, -18}, {62, -14}, {65, -16}, {62, -18}}], 
              Line[{{65, -16}, {67, -16}, {67, 2}, {64, 2}}], 
              Line[{{67, -16}, {67, -19}}]}];
        gr122 = 
          Graphics[{AbsoluteThickness[
                1 + If[i > 1, If[IntegerQ[WCOR[[15]]], WCOR[[15]], 0], 0]], 
              If[i == 1, Line[{{59, -20}, {66, -20}}], 
                Line[{{59, -20}, {61, -20}}]], 
              If[i > 1, Line[{{64, -20}, {66, -20}}], 
                Line[{{59, -20}, {59, -20}}]], 
              If[i > 1, Text[WCOR[[15]], {62.5, -20}], Text["", {0, 0}]]}]; 
        con = WCOR[[15]];
        gr133 = 
          Graphics[{AbsoluteThickness[
                1 + If[i > 1, 
                    Mod[If[IntegerQ[WCOR[[15]]], WCOR[[15]], 0] + crt, 2], 
                    0]], Circle[{67, -20}, 1], Text["+", {67, -20}], 
              Line[{{68, -20}, {70, -20}}], 
              If[i == 1 || i == 3, Line[{{70, -20}, {70, -22}, {71, -23}}], 
                Line[{{70, -20}, {70, -25}, {10, -25}, {10, -21}, {19, 
-21}}]]}];
        gr144 = 
          Graphics[{Table[
                Line[{{4 + 8*(i - 1), -6}, {4 + 8*(i - 1), -8}}], {i, 7}], 
              Line[{{70, -23}, {70, -25}, {10, -25}, {10, -21}, {19, -21}}], 
              Line[{{73, -20}, {75, -20}}]}];
        gr155 = 
          Graphics[{AbsoluteThickness[
                1 + If[i > 1, 
                    Mod[If[IntegerQ[WCOR[[15]]], WCOR[[15]], 0] + crt, 2], 
                    0]], If[i == 3, Line[{{70, -20}, {75, -20}}], 
                Line[{{70, -20}, {72, -20}, {73, -19}}]]}];
        Do[
          If[IntegerQ[WCOR[[i]]], 
            If[WCOR[[i]] == 0, wc = StringInsert[wc, "0", -1], 
              wc = StringInsert[wc, "1", -1]], 
            wc = StringInsert[wc, " ", -1]], {i, 15}];
        If[i > 1, 
          WCOR = ReplacePart[RotateLeft[WCOR, -1], 
              If[i == 2, Mod[WCOR[[15]] + crt, 2], " "], 1]];
        pr = ""; 
        If[i > 1, 
          Do[If[Pnew[[i]] == 0, pr = StringInsert[pr, "0  ", -1], 
              pr = StringInsert[pr, "1  ", -1]], {i, 7}]];
        grtext = 
          Graphics[{Table[Text[Pnew[[i]], {4 + 8*(i - 1), -2}], {i, 8}], 
              Text["&", {63, -16}], Text[">", {30, 7}], Text["<", {60, 2}], 
              Text["<", {10, 2}], Text["v", {0, -5}], Text["<", {30, -25}], 
              Text["^", {67, -5}], Text[">", {40, -17}], 
              Text[If[i == 3, Mod[con + crt, 2], ""], {75, -20}, {-1, 0}], 
              Text["v", {67, -17}], 
              Text[If[i == 1, w[[16 - j]], ""], {-3, 7}, {1, 0}], 
              Text[pr, {30, -10}], 
              Text[If[i == 1, StringTake[rw, -j], 
                  StringTake[wc, 14]], {40, -20}], 
              If[i == 2, Text[Mod[con + crt, 2], {21, -20}], 
                Text["", {0, 0}]]}];
        cr = 
          Graphics[{RGBColor[crt, 0, 0], 
              If[crt == 1, Text["CORRECTION", {30, -12}], 
                Text["", {0, 0}]]}];
        graf[[i, j]] = {gr1, gr2, gr3, gr4, gr5, gr6, gr7, gr8, gr9, gr10, 
            gr11, gr122, gr133, gr144, gr155, grtext, cr}; P = Pnew;
        DWD[[16 - j]] = Mod[con + crt, 2]; crt = correct;
        wc = "", {j, 15}], {i, 3}];
    Do[Do[Print[i, "  ", j];
        Show[graf[[i, j]], PlotRange -> All], {j, 1, 15}], {i, 1, 3}]; 
   If[P == {0, 0, 0, 0, 0, 0, 0, 0}, 
      If[15 == 15 && 3 == 3, 
        Print["THE RECEIVED WORD = ", StringReverse[rw]];
        Do[
          If[DWD[[i]] == 0, dw = StringInsert[dw, "0", -1], 
            dw = StringInsert[dw, "1", -1]], {i, 15}];
        Print["THE DECODE WORD  =  ", StringReverse[dw]]], 
      Print["To decode isn`t possible, we have more than       2errors"]]]

ErrorTrappingDecodeHammingCode[ww_List, x_Symbol] := 
    Module[{cg, k, n, rw, P, gr1, gr1AA, gr2, gr3, vh, gr4, gr5, gr6, gr7, 
        gr8, gr9, gr10, gr11, gr1BB, dcw, DW, w}, w = Reverse[ww]; 
      m = Log[2, Length[w] + 1];
      If[m > 7, Print["The parameter of Hamming code is large"];
        Print["I need the received word of length 256,"];
        Print["not possible to show the process of decoding in the screen."];
        Break[]]; 
      Which[m == 2, g = 1 + x + x^2, m == 3, g = 1 + x + x^3, m == 4, 
        g = 1 + x + x^4, m == 5, g = 1 + x^2 + x^5, m == 6, g = 1 + x + x^6, 
        m == 7, g = 1 + x + x^7];
      If[TEST2[w] == 1, Print["Check the received word on 0 or 1."]; 
        Break[]];
      Print["THE GENERATOR POLYNOMIAL g = ", g, "."];
      cg = CoefficientList[g, {x}]; k = 2^m - m - 1; n = 2^m - 1;
      If[Length[w] != 2^m - 1, 
        Print["The length of the received word = ", Length[w], " ,but  != ", 
          2^m - 1, "."]; Break[]]; rw = "";
      Do[If[w[[i]] == 0, rw = StringInsert[rw, "0", -1], 
          rw = StringInsert[rw, "1", -1]], {i, n}]; P = Table["", {n - k}]; 
      W = P;
      gr1 = Table[Graphics[{}], {n - k}]; gr1AA = gr1;
      gr2 = Table[Graphics[{}], {n - k}]; gr3 = Table[Graphics[{}], {n - k}];
      gr3AA = gr3; vh = Graphics[{}]; gr4 = vh; gr4TE = vh; gr5 = vh; 
      gr6 = vh;
      gr7 = vh; gr8 = vh; gr9 = vh; gr10 = vh; gr11 = vh; gr1BB = vh;
      dcw = Table[" ", {n}]; DW = ""; t = "1";
      Do[t = StringInsert[t, "0", -1], {m - 1}];
      Print["THE RECEIVED WORD = ", StringReverse[rw]];
      Do[Do[Print[f, "  ", j]; 
          gr1BB = Graphics[{AbsoluteThickness[
                  1 + If[IntegerQ[P[[n - k]]], P[[n - k]], 0]], 
                Text["^", {-2 + 8*(n - k), 1.5}], 
                Line[{{-4 + 8*(n - k), 0}, {-2 + 8*(n - k), 
                      0}, {-2 + 8*(n - k), 2}, {-1, 2}, {-1, 1}}], 
                Text["<", {0, 2}], Text[">", {1, 0}]}]; 
          Do[If[cg[[i]] == 1, 
              gr2[[i]] = 
                Graphics[{AbsoluteThickness[
                      1 + If[IntegerQ[P[[n - k]]], P[[n - k]], 0]], 
                    Line[{{-9 + 8*i, 2}, {-9 + 8*i, 1}}], 
                    Text["v", {-9 + 8*i, 1.5}]}]], {i, 2, n - k}];
          
          Do[If[cg[[i]] == 1, 
              If[i == 1, 
                If[f == 1, 
                  W[[i]] = 
                    Mod[w[[n + 1 - j]] + 
                        If[IntegerQ[P[[n - k]]], P[[n - k]], 0], 2], 
                  W[[i]] = If[IntegerQ[P[[n - k]]], P[[n - k]], 0]], 
                W[[i]] = 
                  If[IntegerQ[P[[i - 1]]], 
                    If[IntegerQ[P[[n - k]]], Mod[P[[i - 1]] + P[[n - k]], 2], 
                      P[[i - 1]]], 
                    If[IntegerQ[P[[n - k]]], P[[n - k]], ""]]];
              
              gr1[[i]] = 
                Graphics[{AbsoluteThickness[
                      1 + If[IntegerQ[W[[i]]], W[[i]], 0]], 
                    Text["+", {-9 + 8*i, 0}], Circle[{-9 + 8*i, 0}, 1], 
                    Text[W[[i]], {-5 + 8*i, 0}], 
                    Line[{{-8 + 8*i, 0}, {-6 + 8*i, 0}, {-6 + 8*i, 
                          1}, {-4 + 8*i, 
                          1}, {-4 + 8*i, -1}, {-6 + 8*i, -1}, {-6 + 8*i, 
                          0}}]}], 
              W[[i]] = If[IntegerQ[P[[i - 1]]], P[[i - 1]], ""];
              
              gr1[[i]] = 
                Graphics[{AbsoluteThickness[
                      1 + If[IntegerQ[P[[i - 1]]], P[[i - 1]], 0]], 
                    Text[W[[i]], {-5 + 8*i, 0}], 
                    Line[{{-12 + 8*i, 0}, {-6 + 8*i, 0}, {-6 + 8*i, 
                          1}, {-4 + 8*i, 
                          1}, {-4 + 8*i, -1}, {-6 + 8*i, -1}, {-6 + 8*i, 
                          0}}]}]], {i, 1, n - k}];
          
          Do[gr1AA[[i]] = 
              Graphics[{AbsoluteThickness[1 + If[f == 2, W[[i]], 0]], 
                  Line[{{-5 + 8*i, -1}, {-5 + 8*i, -3}}]}], {i, n - k}];
          
          Do[If[cg[[i + 1]] == 1, 
              gr3[[i]] = 
                Graphics[{AbsoluteThickness[
                      1 + If[IntegerQ[P[[i]]], P[[i]], 0]], 
                    Line[{{-4 + 8*i, 0}, {-2 + 8*i, 0}}]}]], {i, n - k}];
          
          gr4 = Graphics[{AbsoluteThickness[
                  1 + If[f == 1, w[[n + 1 - j]], 0]], 
                Line[{{-10, -10}, {-9, -10}}], 
                Line[{{-2, 0}, {-5, 
                      0}, {-5, -10}, {-6, -10}, {-2, -10}, {-2, -11}, {-2 + 
                        8*(n - k), -11}, {-2 + 
                        8*(n - k), -9}, {-2, -9}, {-2, -10}}]}];
          
          If[f == 1, gr4TE = Graphics[{Text[w[[n + 1 - j]], {-7.5, -10}]}], 
            gr4TE = Graphics[{}]]; 
          gr10 = Graphics[{AbsoluteThickness[
                  1 + If[f == 1, 0, w[[n + 1 - j]]]], 
                Line[{{-2 + 8*(n - k), -10}, {2 + 8*(n - k), -10}}]}];
          
          If[f == 2, 
            If[Table[If[i == 1, 1, 0], {i, n - k}] == W, a1 = 1, a1 = 0]; 
            syn = ""; 
            Do[If[W[[i]] == 0, syn = StringInsert[syn, "0", -1], 
                syn = StringInsert[syn, "1", -1]], {i, n - k}], a1 = 0];
          
          gr7 = Graphics[{AbsoluteThickness[
                  1 + If[f == 2, Mod[a1 + w[[n + 1 - j]], 2], 0]], 
                Circle[{3 + 8*(n - k), -10}, 1], 
                Text["+", {3 + 8*(n - k), -10}], 
                Line[{{4 + 8*(n - k), -10}, {6 + 8*(n - k), -10}}]}];
          
          gr6 = Graphics[{AbsoluteThickness[1 + a1], 
                Line[{{2, -3}, {-2 + 8*(n - k), -3}, {-2 + 
                        8*(n - k), -7}, {2, -7}, {2, -3}}], 
                Line[{{-2 + 8*(n - k), -5}, {3 + 8*(n - k), -5}, {3 + 
                        8*(n - k), -9}}]}];
          
          gr8 = Graphics[{Text["^", {-5, -5}], Text[">", {-3, -10}], 
                Text[">", {8*(n - k), -10}], Text[">", {5 + 8*(n - k), -10}], 
                Text["v", {3 + 8*(n - k), -7}]}];
          
          If[f == 1, 
            If[n - k > 2, 
              gr9 = Graphics[{Text[
                      If[j < 21, StringTake[rw, -j], 
                        StringInsert["...", StringTake[rw, -10], -1]], {-2 + 
                          4*(n - k), -10}], 
                    Text["The circuit of", {4*(n - k), -3.5}], 
                    Text["functional ele-", {4*(n - k), -4.5}], 
                    Text["ments for the", {4*(n - k), -5.5}], 
                    Text["search of errors", {4*(n - k), -6.5}]}], 
              gr9 = Graphics[{Text[
                      If[n <= 12, StringTake[rw, -j], 
                        If[j < 10, StringTake[rw, -j], 
                          StringInsert["...", 
                            StringTake[rw, -9], -1]]], {-2 + 
                          4*(n - k), -10}]}]], 
            If[P == Table[0, {n - k}], 
              gr9 = Graphics[{RGBColor[1, 0, 0], 
                    Text["Syndrome=", {4*(n - k) + 2, -3.5}, {1, 0}], 
                    Text[syn, {4*(n - k) + 2, -3.5}, {-1, 0}], 
                    Text["the received word", {4*(n - k), -4.5}], 
                    Text["is a code word", {4*(n - k), -5.5}]}]; dcw = w;
              vh = Graphics[{Text[w[[n - j + 1]], {8*(n - k) + 8, -10}]}], 
              gr9 = Graphics[{RGBColor[a1, 0, 0], 
                    Text[syn, {4*(n - k) + 2, -3.5}, {-1, 0}], 
                    Text["SYNDROME=", {4*(n - k) + 2, -3.5}, {1, 0}], 
                    Text["is comparing", {4*(n - k), -4.5}], 
                    Text["with", {4*(n - k), -5.5}], 
                    Text[t, {4*(n - k), -6.5}]}];
              
              If[a1 == 1, 
                gr11 = Graphics[{RGBColor[a1, 0, 0], 
                      Text["CORRECTION", {4*(n - k), -8}]}];
                vh = 
                  Graphics[{Text[
                        Mod[w[[n - j + 1]] + 1, 2], {8*(n - k) + 8, -10}]}];
                dcw[[n - j + 1]] = Mod[w[[n - j + 1]] + 1, 2], 
                vh = Graphics[{Text[w[[n - j + 1]], {8*(n - k) + 8, -10}]}];
                dcw[[n - j + 1]] = w[[n - j + 1]]; gr11 = Graphics[{}]]]];
          
          If[f == 2, 
            If[n - k > 2, 
              gr5 = Graphics[{Text[
                      If[Length[w] - j < 22, StringTake[rw, n - j], 
                        StringInsert["...", 
                          StringTake[
                            rw, {n - j - 10, n - j}], -1]], {4*(n - 
                              k), -10}]}], 
              gr5 = Graphics[{Text[
                      If[Length[w] - j < 12, StringTake[rw, n - j], 
                        StringInsert["...", 
                          StringTake[
                            rw, {n - j - 7, n - j}], -1]], {4*(n - k) - 
                          2, -10}]}]]];
          
          Show[vh, gr1, gr1AA, gr2, gr3, gr3AA, gr4, gr4TE, gr5, gr6, gr7, 
            gr8, gr9, gr10, gr11, gr1BB, PlotRange -> All]; 
          P = W, {j, n}], {f, 2}]; 
      Do[If[dcw[[i]] == 0, DW = StringInsert[DW, "0", -1], 
          DW = StringInsert[DW, "1", -1]], {i, n}];
      Print["The received word = ", StringReverse[rw]]; 
      Print["The decode word  =  ", StringReverse[DW]]]

ExpurgatingCode[M_List]:=Prepend[M,Table[1,{Length[M[[1]]]}]]

ExtendedCode[M_List]:=
Prepend[Table[Append[M[[i]],0],{i,1,Length[M]}],
Table[1,{Length[M[[1]]]+1}]]

GeneratorMatrix[H_List]:=NullSpace[H,Modulus->2]

GeneratorPolynomials[n_Integer,x_Symbol] := 
  Module[{FL, FLDL, j, i, l, f, V, gp, g}, If[EvenQ[n],
Union[PolynomialMod[Expand[Flatten[Transpose[
{GeneratorPolynomials[n/2,x]}].{GeneratorPolynomials[n/2,x]}]],2]],
FL = FactorList[x^n - 1, Modulus -> 2]; FLDL = {}; 
     Do[Do[AppendTo[FLDL, FL[[j,1]]], {FL[[j,2]]}], 
      {j, Length[FL]}]; l = Length[FLDL]; 
     If[l > 15, Print["I cannot to find all a generator
 polynomials"]; Break[]]; f = Table[1, {l}]; 
     gp = Table[1, {2^l}]; 
     Do[V = VEC[i, l]; Do[If[V[[j]] == 1, 
         f[[j]] = FLDL[[j]], f[[j]] = 1], {j, l}]; 
       gp[[i + 1]] = 
        PolynomialMod[Expand[Product[f[[i]], {i, l}]], 2] 
       , {i, 2^l - 1}]; g = Union[gp]]]

GilbertVarshamovBound[n_,k_,d_]:=Module[{r,j,i},r=n-k;
If[IntegerQ[n]&&IntegerQ[r]&&IntegerQ[d],
If[Sum[Binomial[n-1,i],{i,0,d-2}]<=2^r,
Print["There exists"];True,
False],
If[IntegerQ[n]&&IntegerQ[d],j=0;
While[Sum[Binomial[n-1,i],{i,0,d-2}]<=2^(n-j),
j=j+1];k<=j-1,
If[IntegerQ[k]&&IntegerQ[d],j=k+d-1;
While[Sum[Binomial[j-1,i],{i,0,d-2}]>2^(j-k),
j=j+1];n>j-1,
If[IntegerQ[n]&&IntegerQ[k],j=2;
While[Sum[Binomial[n-1,i],{i,0,j-2}]<=2^(n-k),
j=j+1];d<=j-1,
Print["If 1 + (n-1)/1! + n(n-1)/2! + ... +
Binomial[n-1,d-2] < ",2^r];
Print["there exists a binary linear code  of length ", n];
Print["with at most ",n-k," parity checks"];
Print["and code  distance at least ",d," ."]]]]]]

HammingBound[n_,k_,d_]:=Module[{s,i},
If[EvenQ[d],Print[" d Odd number in Hamming bound"];Break[]];
If[IntegerQ[n]&&IntegerQ[k]&&IntegerQ[d],
If[(2^k)*Sum[Binomial[n,i],{i,0,Floor[(d-1)/2]}]<=2^n,
True,False],
 If[IntegerQ[n]&&IntegerQ[k],
j=2;While[Sum[(2^k)*Binomial[n,i],{i,0,Floor[(j)/2]}]<=2^n,
j=j+1];d<=If[Sum[(2^k)*Binomial[n,i],{i,0,Floor[(j-1)/2]}]==2^n,j-1,j],
  If[IntegerQ[n]&&IntegerQ[d],
k<=Floor[N[Log[2,(2^n)/Sum[Binomial[n,i],{i,0,Floor[(d)/2]}]]]]
,If[IntegerQ[k]&&IntegerQ[d],s=d+k-1;
While[(2^k)*(Sum[Binomial[s,i],{i,0,Floor[d/2]}])>2^s,
s=s+1];n>=s ,
Print["For t",If[IntegerQ[d]," = ",""],
If[IntegerQ[d],Floor[d/2],""],
" - error - correcting binary code of length ", n];
Print["Hamming bound :  ",
If[IntegerQ[d]&&(d<=9),
(2^k)Sum[Binomial[n,i],{i,0,Floor[d/2]}]<=2^n,
2^k],
If[Not[IntegerQ[d]]||(d>=11),"(1 + n + n(n-1)/2 + ...
+n(n-1)(n-2)...(n-t)/t!) <= ",""],
If[Not[IntegerQ[d]]||(d>=11),2^n,""]
];If[IntegerQ[d]&&(d<=9),(2^k)Sum[Binomial[n,i],{i,0,Floor[d/2]}]<=2^n]]]]]]

HammingCode[m_Integer, inf_List] := 
  Module[{M, i, NS, W}, If[TEST2[inf] == 1, 
     Print["Check of the inf on 0 or 1"]; Break[]]; 
    If[Length[inf] != 2^m - 1 - m, 
     Print["Length of inf = ", Length[inf], "  =!=  ", 
       2^m - 1 - m]; Break[]]; 
    M = HammingMatrix[m]; 
    NS = NullSpace[M, Modulus -> 2]; 
     W = Mod[inf . NS, 2]] 

HammingMatrix[m_Integer] := 
  Transpose[Table[VEC[i, m], {i, 2^m - 1}]]

McWilliamsIdentity[f_,n_Integer,x_Symbol]:=
Expand[Together[((1+x)^n/(f/.x->1))*(f/.x->(1-x)/(1+x))]]

MeggittDecode[g_, ww_List, n_Integer, x_Symbol] := 
  Module[{cg, g1, k, b1, gr1, gr2, gr3, gr4, gr6, gr7, gr8, gr9, gr10, gr11, 
      gr1AA, gr1BB, gr4TE, vh, dcw, DW, inf, d, d1, j1, P, W, f, syn, a1, 
      Dsure, i, j, l, dkon},w=Reverse[ww]; g1 = PolynomialMod[g, 2]; b1 = 0; b2 = 0; 
    dkon = 0;
    l = Length[
        CoefficientList[PolynomialMod[x^n - 1, g1, Modulus -> 2], {x}]]; 
    If[l != 0, 
      Print["The polynomial g =", g1, " isn`t the  generator polynomial,as"];
      Print["The generator polynomial of cyclic code of  length ", n, 
        " must divides ", x^n - 1]; Break[]];
    If[TEST2[w] == 1, Print["Check the received word on 0 or 1."]; Break[]];
    cg = CoefficientList[g1, {x}]; k = n - Length[cg] + 1;
    If[Length[w] != n, 
      Print["The length of the received word = ", Length[w], " ,but  != ", n, 
        "."]; Break[]]; rw = "";
    Do[If[w[[i]] == 0, rw = StringInsert[rw, "0", -1], 
        rw = StringInsert[rw, "1", -1]], {i, n}]; P = Table["", {n - k}]; 
    W = P;
    gr1 = Table[Graphics[PostScript[""]], {n - k}]; gr1AA = gr1;
    gr2 = Table[Graphics[PostScript[""]], {n - k}]; 
    gr3 = Table[Graphics[PostScript[""]], {n - k}];
    gr3AA = gr3; vh = Graphics[PostScript[""]]; gr4 = vh; gr4TE = vh; 
    gr5 = vh; gr6 = vh;
    gr7 = vh; gr8 = vh; gr9 = vh; gr10 = vh; gr11 = vh; gr1BB = vh;
    dcw = Table[" ", {n}]; DW = "";
    If[IntegerQ[Stop], SL = SYNDROMLIST[g, n, k, Distance, x], 
      If[IntegerQ[Dsure], Print["The code distance = ", Distance];
        SL = SYNDROMLIST[g, n, k, Distance, x], 
        If[k < 13, d = n; Do[inf = VEC[j1, k];
            cw = PolynomialMod[Sum[inf[[i]]*x^(i - 1)*g, {i, k}], 2];
            d1 = cw /. x -> 1; d = If[d > d1, d1, d];
            
            If[Mod[j1, 2^8] == 0, 
              Print[j1, " Please wait, I find the code distance"]], {j1, 
              2^k - 1}]; Print["The code distance = ", d, "."]; dkon = d;
          SL = SYNDROMLIST[g, n, k, d, x], 
          If[IntegerQ[Distance], 
            If[Distance <= 2, Print["This code cannot correct"];
              Print["even a single error"]; Clear[Distance, Stop, Dsure];
              Break[]]; 
            If[TEST3[g, n, Distance] == 1, 
              Print["if you want, you can try again"];
              Print["Distance= ;MeggittDecode[g,w,n,x]"]; Break[]];
            If[IntegerQ[Stop], Break[]]; 
            SL = SYNDROMLIST[g, n, k, Distance, x], 
            Print["The efficiency of code is large,"];
            Print["you must give the code distance and try again."];
            Print["Distance=   ;MeggittDecode[g,w,n,x]"]; Break[]]]];
      If[IntegerQ[Stop], Break[]]]; Clear[Stop];
    Print["THE SYNDROME LIST =", MatrixForm[SL]];
    Do[Do[Print[f, "  ", j]; 
        gr1BB = Graphics[{AbsoluteThickness[
                1 + If[IntegerQ[P[[n - k]]], P[[n - k]], 0]], 
              Text["^", {-2 + 8*(n - k), 1.5}], 
              Line[{{-4 + 8*(n - k), 0}, {-2 + 8*(n - k), 0}, {-2 + 8*(n - k),
                     2}, {-1, 2}, {-1, 1}}], Text["<", {0, 2}], 
              Text[">", {1, 0}]}]; 
        Do[If[cg[[i]] == 1, 
            gr2[[i]] = 
              Graphics[{AbsoluteThickness[
                    1 + If[IntegerQ[P[[n - k]]], P[[n - k]], 0]], 
                  Line[{{-9 + 8*i, 2}, {-9 + 8*i, 1}}], 
                  Text["v", {-9 + 8*i, 1.5}]}]], {i, 2, n - k}];
        Do[
          If[cg[[i]] == 1, 
            If[i == 1, 
              If[f == 1, 
                W[[i]] = 
                  Mod[w[[n + 1 - j]] + 
                      If[IntegerQ[P[[n - k]]], P[[n - k]], 0], 2], 
                W[[i]] = If[IntegerQ[P[[n - k]]], P[[n - k]], 0]], 
              W[[i]] = 
                If[IntegerQ[P[[i - 1]]], 
                  If[IntegerQ[P[[n - k]]], Mod[P[[i - 1]] + P[[n - k]], 2], 
                    P[[i - 1]]], If[IntegerQ[P[[n - k]]], P[[n - k]], ""]]];
            
            gr1[[i]] = 
              Graphics[{AbsoluteThickness[
                    1 + If[IntegerQ[W[[i]]], W[[i]], 0]], 
                  Text["+", {-9 + 8*i, 0}], Circle[{-9 + 8*i, 0}, 1], 
                  Text[W[[i]], {-5 + 8*i, 0}], 
                  Line[{{-8 + 8*i, 0}, {-6 + 8*i, 0}, {-6 + 8*i, 
                        1}, {-4 + 8*i, 
                        1}, {-4 + 8*i, -1}, {-6 + 8*i, -1}, {-6 + 8*i, 0}}]}],
             W[[i]] = If[IntegerQ[P[[i - 1]]], P[[i - 1]], ""];
            
            gr1[[i]] = 
              Graphics[{AbsoluteThickness[
                    1 + If[IntegerQ[P[[i - 1]]], P[[i - 1]], 0]], 
                  Text["v", {-9 + 8*i, -2}], Text[W[[i]], {-5 + 8*i, 0}], 
                  Line[{{-12 + 8*i, 0}, {-6 + 8*i, 0}, {-6 + 8*i, 
                        1}, {-4 + 8*i, 
                        1}, {-4 + 8*i, -1}, {-6 + 8*i, -1}, {-6 + 8*i, 
                        0}}]}];
            
            gr1AA[[i]] = 
              Graphics[{AbsoluteThickness[
                    1 + If[f == 2, If[IntegerQ[P[[i - 1]]], P[[i - 1]], 0], 
                        0]], Line[{{-9 + 8*i, 0}, {-9 + 8*i, -3}}]}]], {i, 1, 
            n - k}];
        Do[
          If[cg[[i + 1]] == 1, 
            gr3[[i]] = 
              Graphics[{AbsoluteThickness[
                    1 + If[IntegerQ[P[[i]]], P[[i]], 0]], 
                  Line[{{-4 + 8*i, 0}, {-2 + 8*i, 0}}]}];
            
            gr3AA[[i]] = 
              Graphics[{AbsoluteThickness[
                    1 + If[f == 2, If[IntegerQ[P[[i]]], P[[i]], 0], 0]], 
                  Line[{{-3 + 8*i, 0}, {-3 + 8*i, -3}}], 
                  Text["v", {-3 + 8*i, -2}]}]], {i, n - k}];
        gr4 = 
          Graphics[{AbsoluteThickness[1 + If[f == 1, w[[n + 1 - j]], 0]], 
              Line[{{-10, -10}, {-9, -10}}], 
              Line[{{-2, 0}, {-5, 
                    0}, {-5, -10}, {-6, -10}, {-2, -10}, {-2, -11}, {-2 + 
                      8*(n - k), -11}, {-2 + 
                      8*(n - k), -9}, {-2, -9}, {-2, -10}}]}];
        If[f == 1, gr4TE = Graphics[{Text[w[[n + 1 - j]], {-7.5, -10}]}], 
          gr4TE = Graphics[PostScript[""]]]; 
        gr10 = Graphics[{AbsoluteThickness[1 + If[f == 1, 0, w[[n + 1 - j]]]],
               Line[{{-2 + 8*(n - k), -10}, {2 + 8*(n - k), -10}}]}];
        If[f == 2, 
          Do[If[SL[[i]] == P, a1 = 1; b1 = b1 + 1; Break[], a1 = 0], {i, 
              Length[SL]}];
          
          If[b1 > t, Print["Excuse me,the code distance less than ", dis]; 
            Print["you can try again"];
            Print["Distance=  ;MeggittDecode[g,w,n,x]"]; Break[]];
          syn = ""; 
          Do[If[P[[i]] == 0, syn = StringInsert[syn, "0", -1], 
              syn = StringInsert[syn, "1", -1]], {i, n - k}], a1 = 0];
        gr7 = 
          Graphics[{AbsoluteThickness[
                1 + If[f == 2, Mod[a1 + w[[n + 1 - j]], 2], 0]], 
              Circle[{3 + 8*(n - k), -10}, 1], 
              Text["+", {3 + 8*(n - k), -10}], 
              Line[{{4 + 8*(n - k), -10}, {6 + 8*(n - k), -10}}]}];
        gr6 = 
          Graphics[{AbsoluteThickness[1 + a1], 
              Line[{{2, -3}, {-2 + 8*(n - k), -3}, {-2 + 
                      8*(n - k), -7}, {2, -7}, {2, -3}}], 
              Line[{{-2 + 8*(n - k), -5}, {3 + 8*(n - k), -5}, {3 + 
                      8*(n - k), -9}}]}];
        gr8 = 
          Graphics[{Text["^", {-5, -5}], Text[">", {-3, -10}], 
              Text[">", {8*(n - k), -10}], Text[">", {5 + 8*(n - k), -10}], 
              Text["v", {3 + 8*(n - k), -7}]}];
        If[f == 1, 
          If[n - k > 2, 
            gr9 = Graphics[{Text[
                    If[j < 21, StringTake[rw, -j], 
                      StringInsert["...", StringTake[rw, -10], -1]], {-2 + 
                        4*(n - k), -10}], 
                  Text["The circuit of", {4*(n - k), -3.5}], 
                  Text["functional ele-", {4*(n - k), -4.5}], 
                  Text["ments for the", {4*(n - k), -5.5}], 
                  Text["search of errors", {4*(n - k), -6.5}]}], 
            gr9 = Graphics[{Text[
                    If[n <= 12, StringTake[rw, -j], 
                      If[j < 10, StringTake[rw, -j], 
                        StringInsert["...", StringTake[rw, -9], -1]]], {-2 + 
                        4*(n - k), -10}]}]], 
          If[P == Table[0, {n - k}], b2 = 1; 
            gr9 = Graphics[{RGBColor[1, 0, 0], 
                  Text["Syndrome=", {4*(n - k), -3.5}, {1, 0}], 
                  Text[syn, {4*(n - k), -3.5}, {-1, 0}], 
                  Text["the received word", {4*(n - k), -4.5}], 
                  Text["is a code word", {4*(n - k), -5.5}]}]; dcw = w;
            vh = Graphics[{Text[w[[n - j + 1]], {8*(n - k) + 8, -10}]}], 
            gr9 = Graphics[{Text[syn, {4*(n - k), -3.5}], 
                  Text["SYNDROME", {4*(n - k), -4.5}]}];
            
            If[a1 == 1, 
              gr11 = Graphics[{RGBColor[a1, 0, 0], 
                    Text["CORRECTION", {4*(n - k), -5.5}]}];
              
              vh = Graphics[{Text[
                      Mod[w[[n - j + 1]] + 1, 2], {8*(n - k) + 8, -10}]}];
              dcw[[n - j + 1]] = Mod[w[[n - j + 1]] + 1, 2], 
              vh = Graphics[{Text[w[[n - j + 1]], {8*(n - k) + 8, -10}]}];
              dcw[[n - j + 1]] = w[[n - j + 1]]; 
              gr11 = Graphics[PostScript[""]]]]];
        If[f == 2, 
          If[n - k > 2, 
            gr5 = Graphics[{Text[
                    If[Length[w] - j < 22, StringTake[rw, n - j], 
                      StringInsert["...", 
                        StringTake[
                          rw, {n - j - 10, n - j}], -1]], {4*(n - k), -10}]}],
             gr5 = Graphics[{Text[
                    If[Length[w] - j < 12, StringTake[rw, n - j], 
                      StringInsert["...", 
                        StringTake[
                          rw, {n - j - 7, n - j}], -1]], {4*(n - k) - 
                        2, -10}]}]]];
        Show[vh, gr1, gr1AA, gr2, gr3, gr3AA, gr4, gr4TE, gr5, gr6, gr7, gr8, 
          gr9, gr10, gr11, gr1BB, PlotRange -> All]; P = W, {j, n}], {f, 2}]; 
    Do[If[dcw[[i]] == 0, DW = StringInsert[DW, "0", -1], 
        DW = StringInsert[DW, "1", -1]], {i, n}]; If[b1 > t, Break[]];
    If[b2 == 1 || b1 > 0, Print["THE RECEIVED WORD = ", StringReverse[rw]];
      Print["THE DECODED WORD  = ", StringReverse[DW]], 
      If[dkon == 0, Print["This received word cannot decode "];
        Print["either you give a wrong code distance ,"];
        Print["or we have more than ", t, " error", If[t > 1, "s", ""]], 
        Print["This received word cannot decode ,"];
        Print["as we have more than ", t, " error", If[t > 1, "s", ""]]]]]

NonsystematicEncodeCyclicCode[g_, ww_List, n_Integer, x_Symbol] := 
  Module[{g1, cg, u, v, M, L, k1, d, i, j, cdw, k, a, b, T, W, P, ph, gr, A, 
      B, inf},
    inf = Reverse[ww]; 
    If[TEST2[inf] == 1, Print["Check inf on 0 or 1"]; Break[]];
    g1 = PolynomialMod[g, 2]; u = Sum[inf[[i]]*x^(i - 1), {i, Length[inf]}];
    l = Length[
        CoefficientList[PolynomialMod[x^n - 1, g1, Modulus -> 2], {x}]]; 
    If[l != 0, 
      Print["The polynomial g= ", g1, " isn`t a generator polynomial, as"];
      Print["the generator polynomial of cyclic code of length ", n, 
        " must divides ", x^n - 1]; Break[]];
    cg = Reverse[CoefficientList[g1, {x}]]; k1 = n - Length[cg] + 1;
    v = PolynomialMod[Expand[g1*u], 2]; M = Reverse[CoefficientList[v, {x}]];
    If[Length[M] < n, Do[PrependTo[M, 0], {n - Length[M]}]];
    L = Reverse[CoefficientList[u, {x}]];
    If[Length[inf] != k1, 
      Print["The length of inf =", Length[inf], ",but must to be = ", k1]; 
      Break[]];
    If[Length[L] < k1, Do[PrependTo[L, 0], {k1 - Length[L]}]];
    Do[L = Insert[L, 0, -1], {n - k1}]; d = ""; cdw = "";
    Do[If[L[[i]] == 0, d = StringInsert[d, "0", -1], 
        d = StringInsert[d, "1", -1]], {i, n}]; d = StringTake[d, k1];
    Do[d = StringInsert[d, " ", -1], {n - k1}];
    Print["INFORMATION POLYNOMIAL=", u];
    Print["INFORMATION WORD =", 
      StringReverse[StringInsert[StringReverse[d], " ", n - k1 + 1]]];
    Do[If[M[[i]] == 0, cdw = StringInsert[cdw, "0", -1], 
        cdw = StringInsert[cdw, "1", -1]], {i, n}];
    cdw = StringInsert[cdw, " ", 1]; k = Length[cg]; P = Table[0, {k - 1}];
    P1 = Table[" ", {k - 1}]; gr = Table[Graphics[{}], {k}]; ph = Table[Graphics[{}], {k - 1}];
    a = Table[0, {n}, {k}]; b = Table[0, {n}, {k}];
    Do[W = P; Do[W[[k - i]] = W[[k - i - 1]], {i, k - 2}]; W[[1]] = L[[j]];
      a[[j, 1]] = L[[j]]; Do[a[[j, i]] = P[[i - 1]], {i, 2, k}]; T = P; 
      P = W;
      If[j < k - 1, 
        Do[P1 = Flatten[ReplacePart[P1, Take[P, {i, i}], i]], {i, j}], 
        P1 = P]; 
      If[j > k1, P2 = P; 
        Do[P1 = Flatten[ReplacePart[P2, " ", i]]; P2 = P1; 
          Null, {i, j - k1}]]; 
      gr[[1]] = 
        Graphics[{AbsoluteThickness[1 + a[[j, 1]]], 
            Line[{{-4 - 1.1*(n - j), 0}, {-5, 0}}], 
            Line[{{-1, 0}, {0, 0}, {0, 5}, {3, 5}}], Line[{{0, 0}, {1, 0}}], 
            Line[{{1, 1}, {3, 1}, {3, -1}, {1, -1}, {1, 1}}], 
            Text[P1[[1]], {2, 0}]}];
      Do[gr[[i]] = 
          If[cg[[i]] == 0, 
            Graphics[{AbsoluteThickness[1 + a[[j, i]]], 
                Line[{{4*i - 5, 0}, {4*i - 3, 0}}], 
                Line[{{4*i - 3, 1}, {4*i - 1, 
                      1}, {4*i - 1, -1}, {4*i - 3, -1}, {4*i - 3, 1}}], 
                Text[P1[[i]], {4*i - 2, 0}]}], 
            Graphics[{AbsoluteThickness[1 + a[[j, i]]], 
                Line[{{4*i - 5, 0}, {4*i - 3, 0}}], 
                Line[{{4*i - 3, 1}, {4*i - 1, 
                      1}, {4*i - 1, -1}, {4*i - 3, -1}, {4*i - 3, 1}}], 
                Line[{{4*i - 4, 0}, {4*i - 4, 4}}], 
                Text[P1[[i]], {4*i - 2, 0}]}]], {i, 2, k - 1}];
      gr[[k]] = 
        Graphics[{AbsoluteThickness[1 + a[[j, k]]], 
            Line[{{4*(k - 1) - 1, 0}, {4*(k - 1), 0}, {4*(k - 1), 4}}]}];
      Do[ph[[i]] = If[cg[[i + 1]] == 0, b[[j, 1]] = a[[j, 1]];
            b[[j, i + 1]] = b[[j, i]];
            
            Graphics[{AbsoluteThickness[1 + b[[j, i + 1]]], 
                Line[{{4*i - 1, 5}, {4*i + 3, 5}}]}], b[[j, 1]] = a[[j, 1]]; 
            b[[j, i + 1]] = Mod[b[[j, i]] + T[[i]], 2];
            
            Graphics[{AbsoluteThickness[1 + b[[j, i + 1]]], 
                Circle[{4*i, 5}, 1], Text["+", {4*i, 5}], 
                Line[{{4*i + 1, 5}, {4*i + 3, 5}}]}]], {i, k - 2}];
      If[k > 2, b[[j, k]] = Mod[b[[j, k - 1]] + a[[j, k]], 2], 
        b[[j, k]] = Mod[a[[j, 1]] + T[[1]], 2]];
      ph[[k - 1]] = 
        Graphics[{AbsoluteThickness[1 + b[[j, k]]], Circle[{4*(k - 1), 5}, 1],
             Text["+", {4*(k - 1), 5}], Line[{{4*k - 3, 5}, {4*k - 1, 5}}], 
            Line[{{4*k + 3, 5}, {4*k + 1.1*j + 2, 5}}]}];
      A = 
        Graphics[{Text[StringTake[d, {j}], {-2, 0}, {1, 0}], 
            Text[StringTake[StringReverse[d], n - j], {-5, 0.5}, {1, 0}]}];
      B = 
        Graphics[{Text[StringTake[cdw, {j + 1}], {4*k, 5}, {-1, 0}], 
            Text[StringReverse[StringTake[cdw, j]], {4*k + 3, 5.5}, {-1, 
                0}]}];
      Show[ph, gr, A, B], {j, n}]; Print["CODE POLYNOMIAL=", v];
Print["CODE WORD = ", cdw]]

ParityCheckMatrix[G_List]:=NullSpace[G,Modulus->2]

PuncturedCode[M_List,n_Integer]:=
ParityCheckMatrix[
Transpose[Drop[Transpose[GeneratorMatrix[M]],{n,n}]]]

RectangularCode[n_Integer, m_Integer, inf_List] := 
  Module[{a, i, j, w, b, c}, 
   If[TEST2[inf] == 1, Print["CHECK inf on 0 or 1 "]; 
      Break[]]; If[Length[inf] != n*m, 
     Print["CHECK THE LENGTH inf"]; Break[]]; a = inf; 
    Do[a = Insert[a, X, i*(n + 1)], {i, m}]; 
    Do[AppendTo[a, X], {i, n + 1}]; a = Partition[a, n + 1];
    Do[a[[i,n + 1]] = sum2[Drop[a[[i]], {n + 1}]], {i, m}]; 
    Do[a[[m + 1,j]] = sum2[Table[a[[i,j]], {i, m}]], 
     {j, n + 1}]; w = Flatten[a]]

ShowBinaryGaloisField[g_,x_Symbol,b_Symbol,a_Symbol] := 
  Module[{n, o, od, br, pr, W, KL, j, i, prv, r, 
    f, v, len, minpol, probel, FPOL,Pole},Pole={};
   n = Exponent[g, x]; 
    If[n == 1, Print["Cive the another\
        polynomial with degree more then 1"]; 
      Break[]]; FL = FactorList[g, Modulus -> 2]; 
    If[Length[FL] != 2 || FL[[2,2]] != 1, 
     Print["Polynomial ", g, 
       " is not irreducible polynomial"]; Break[]
      ];If[n>=10,Print["The Galois Field has ",2^n," elements"];
Print["Use the program"];
Print["BinaryGaloisField[",InputForm[g],",x,b,a]"];Break[]];  o = ""; Do[o = 
      StringInsert[o, "0", -1], {n}]; 
    od = StringTake[StringInsert[o, "1", 1], n]; 
    br = 0; pr = " "; W = 0; 
    Do[If[br == 0, v = Reverse[VEC[i, n]]; 
       f = Sum[v[[j]]*x^(j - 1), {j, n}]; 
       Do[If[Mod[2^n - 1, j] == 0, 
         If[j != 2^n - 1, 
           r = 
            PolynomialMod[f^j, g, 
             Modulus -> 2]]; 
          If[r == 1 && j < 2^n - 1, Break[], 
           pr = f; br = 1]], {j, 2^n - 1}]], 
     {i, 2^n - 1}]; Do[f = 
       PolynomialMod[pr^i, g, Modulus -> 2] /. 
        x -> a; FPOL = 
       FactorList[x^(2^n - 1) - 1, Modulus -> 2]\
       ; stop = 0; Do[If[stop == 1, Break[]]; 
        ar = 
         AlgebraicRules[(g /. x -> a) == 0, 
          {x, a}]; 
        If[PolynomialMod[PolynomialMod[FP\
               OL[[j,1]] /. x -> f, 2] /. ar\
            , 2] == 0, 
         minpol = FPOL[[j,1]]; stop = 1; 
          Break[]], {j, Length[FPOL]}]; 
      KL = CoefficientList[f, {a}]; prv = ""; 
      Do[If[KL[[i]] == 0, 
        prv = StringInsert[prv, "0", -1], 
        prv = StringInsert[prv, "1", -1]], 
       {i, Length[KL]}]; 
      len = StringLength[prv]; 
      If[len < n, Do[prv = 
         StringInsert[prv, "0", -1], {n - len}]]\
       ; If[W == 0, 
       Print[" "]; Print["    Galois Field GF(", 2^n, 
         ")"];Print[" with irreducible  polynopmial "];
 Print["              ",g];
Print["____________________________________________"];
 probel = "  "; 
        Do[probel = 
          StringInsert[probel, " ", 1], 
         {5*(n - 1) - POLPRINTLENGTH[1,x]}];
Print["Log ","Vector", " ", 
         "  pr.el.","polynom" , StringTake[probel,StringLength[probel]-1], "min.polynom"]; 
Print["____________________________________________"]; 
        Print["-",\[Infinity], "    ", o, "   ", 0, 
         "      ", 0, probel, x]; 
AppendTo[Pole,{-\[Infinity], o, 0, 0, x}];
        Print["  ", 0, "    ", od, "   ", 1, 
         "      ", 1, probel, 1 + x];
AppendTo[Pole,{0,od,1,1,1 + x}]];
 W = 1; 
      probel = "  ";
      Do[probel = StringInsert[probel, " ", 1], 
       {5*n - 5 - POLPRINTLENGTH[f /. a -> x,x]}]; 
      Print["  ", i, If[i < 10, "    ", "  "], 
       prv, "   ", b^i, 
       Which[i == 1, "      ", i > 1 && i < 10, 
        "     ", i > 9 && i < 100, "    ", 
        i > 99, "  "], f, probel, minpol];
AppendTo[Pole,{i,prv,b^i,f,minpol}], 
     {i, 1, 2^n - 2}]]

ShowBurstRectangularCodeVector[n_Integer, m_Integer, inf_List,X_Symbol] := 
  Module[{a, i, j, w, b, c}, If[TEST2[inf] == 1, 
     Print["CHECK VECTOR inf ON 0 OR 1"]; Break[]]; a = inf; 
    If[Length[a] != n*m, Print["THE LENGTH INF =!= n*m"]; Break[]]; 
    Print["INFORM=", a]; Print["  "]; Do[a = Insert[a, X, i*(n + 1)], {i, m}]; 
    Do[AppendTo[a, X], {i, n + 1}]; a = Partition[a, n + 1]; 
    Print["MATRIX FORM=", MatrixForm[a]]; Print["  "]; 
    Do[a[[i,n + 1]] = sum2[Drop[a[[i]], {n + 1}]], {i, m}]; 
    Do[a[[m + 1,j]] = sum2[Table[a[[i,j]], {i, m}]], {j, n + 1}]; 
    Print["MATRIX FORM=", MatrixForm[a]]; Print[" "]; 
    w = Table[0, {(n + 1)*(m + 1)}]; 
    Do[Do[w[[Mod[(m + 1)*(j - i) + i - 1, (m + 1)*(n + 1)] + 1]] = a[[i,j]], 
      {i, m + 1}], {j, n + 1}]; Print["CODE WORD=", w]]

ShowCorrectBurstRectangularCode[n_Integer, m_Integer, w_List,X_Symbol] := 
  Module[{y, A, B, BT, h, v, er, i, j, k, u, i0}, 
   If[Length[w] != (n + 1)*(m + 1), 
     Print["CHECK THE LENGTH w=", Length[w], "!=", (n + 1)*(m + 1)]; Break[]]; 
    If[TEST2[w] == 1, Print["CHECK w ON 0 OR 1"]; Break[]]; 
    If[n + 1 < 2*(m + 1) - 3, Print["THE RECTANGULAR CODE IS"];
Print["m-BURST CORRECTING IF AND ONLY IF n+1>=2*(m+1)-3"]; Break[]]; Print["RECEIVED WORD=", w]; 
    Print[" "]; A = Table[0, {m + 1}, {n + 1}]; 
    Do[Do[A[[i,j]] = w[[Mod[(m + 1)*(j - i) + i - 1, (m + 1)*(n + 1)] + 1]], 
      {i, m + 1}], {j, n + 1}]; Print[" "]; h = Table[0, {m + 1}]; 
    v = Table[0, {n + 1}]; Do[er = sum2[A[[i]]]; 
      If[er == 0, Continue[], h[[i]] = 1], {i, m + 1}]; 
    Print["MATRIX FORM=", MatrixForm[A], "   SYNDROME h  ", TableForm[h]]; 
    Do[er = sum2[Table[A[[i,j]], {i, m + 1}]]; 
      If[er == 0, Continue[], v[[j]] = 1], {j, n + 1}]; Print[" "]; 
    Print[" SYNDROME v=", v]; If[Apply[Plus, h] == 0 && Apply[Plus, v] == 0, 
     Print["THE RECEIVED WORD IS A CODE WORD"]; Break[]]; 
    If[Apply[Plus, h] != Apply[Plus, v], 
     Print["DECODING ISN`T POSSIBLE"]; Break[]]; 
    u = Table[{i, 0, 0}, {i, n + 1}]; 
    Do[If[v[[i]] == 1 && v[[PP[i + 1, n + 1]]] == 0, 
      i0 = i; k = i + 1; j = 0; 
       While[v[[PP[k, n + 1]]] == 0, j += 1; k += 1]; 
       u[[i]] = {i0, j, PP[k, n + 1]}], {i, n + 1}]; 
    For[i = 1, u[[i]][[2]] < n - m + 1 && i <= n + 1, i++, 
     If[i > n, i = i + 1; Break[]]]; 
    If[i > n + 1, Print["DECODING ISN`T POSSIBLE"]; Break[]]; y = u[[i]]; 
    B = Table[0, {m + 1}, {n + 1}]; j = y[[3]]; 
    Do[If[h[[i]] == 1, B[[i,PP[j, n + 1]]] = 1; Label[START]; j++; 
       If[v[[PP[j, n + 1]]] == 0, Goto[START]]], {i, m + 1}]; Print[" "]; 
    Print["THE BURST =", MatrixForm[B]]; Print[" "]; 
    Print["  THE RECEIVED WORD=", w]; BT = Table[0, {(n + 1)*(m + 1)}]; 
    Do[Do[BT[[Mod[(m + 1)*(j - i) + i - 1, (m + 1)*(n + 1)] + 1]] = 
       B[[i,j]], {i, m + 1}], {j, n + 1}]; 
    If[Apply[Max, Flatten[Position[BT, 1]]] - 
       Apply[Min, Flatten[Position[BT, 1]]] > m, 
     Print["DECODING ISN`T POSSIBLE"]; Break[]]; 
    Print["      THE BURST  =  ", BT]; 
    Print["THE CORRECTING WORD=", Mod[BT + w, 2]]]

ShowCyclicCode[g_, n_Integer,x_Symbol] := 
  Module[{g1, y, Y, l, k, K, G, i, j, R, H},   
    g1 = PolynomialMod[g, 2]; y = 0; 
    Y = 0; If[n > 10, y = 1; Y = 2]; 
    If[n > 76, Print["The length of code is large,"]; 
      Print["the generator and parity check matrix"]; 
      Print["will not be shown on the screen."]; 
      Print["You can use command:"];
Print["CyclicCode[",g,",",n,"]"];Break[]]; 
    l = Length[CoefficientList[PolynomialMod[x^n - 1, g1, 
        Modulus -> 2], {x}]]; 
    If[l != 0, Print["The polynomial g =", g1, 
       " isn`t a generator polynomial, as"]; 
      Print["the generator polynomial of cyclic code of
length ", n, " must divides  ", x^n - 1]; Break[]]; 
    Print["   GENERATOR"]; 
    Print["   POLYNOMIAL     g = ", g1]; Print[" "]; 
    k = n - Length[CoefficientList[g1, {x}]] + 1; 
    GM = Array[a, {k, n}]; 
    Do[L = CoefficientList[g1*x^(i - 1), {x}]; 
      Do[AppendTo[L, 0], {k - i}]; GM[[i]] = L, {i, k}]; 
    If[y == 0, Print["GENERATOR MATRIX GM = ", 
       MatrixForm[GM]]; Print[" "]]; 
    If[Y == 2, K = 
       {"0", "1"}; 
      G = Table["", {k}]; 
      Do[Do[G[[i]] = 
         StringInsert[G[[i]], K[[GM[[i,j]] + 1]], -1], 
        {j, n}], {i, k}]; 
      Print["GENERATOR MATRIX GM =", MatrixForm[G]]]; 
    h = PolynomialMod[PolynomialQuotient[x^n - 1, g1, x], 
      2]; Print["PARITY"]*
     Print["CHECK POLYNOMIAL h =  ", h]; 
    HM = Array[b, {n - k, n}]; 
    Do[R = Reverse[CoefficientList[h, {x}]]; 
      Do[AppendTo[R, 0], {n - k - i}]; 
      Do[PrependTo[R, 0], {i - 1}]; HM[[n - k - i + 1]] = R\
      , {i, n - k}]; If[y == 0, 
     Print[" "]; Print["PARITY "]; 
      Print["CHECK MATRIX   HM =  ", MatrixForm[HM]]]; 
    If[Y == 2, H = Table["", {n - k}]; 
      Do[Do[H[[i]] = 
         StringInsert[H[[i]], K[[HM[[i,j]] + 1]], -1], 
        {j, n}], {i, n - k}]; Print[" "]; Print["PARITY"]; 
      Print["CHECK MATRIX   HM =  ", MatrixForm[H]]]]

ShowDecHammingCode[m_Integer, w_List] := 
  Module[{M, NS, ER, S}, 
   If[TEST2[w] == 1, Print["Check w on 0 or 1"]; Break[]]; 
    If[Length[w] != 2^m - 1, 
     Print["Length w = ", Length[w], " =!= ", 2^m - 1]; 
      Break[]]; If[m > 4,Print["The parameter m is large"]; 
      Print["The parity check matrix will not be show on
the screen"]]; M = HammingMatrix[m]; 
    S = Mod[M . w, 2]; If[m <= 4, 
     Print["PARITY"]; Print["CHECK"]; 
      Print["MATRIX=", MatrixForm[M], "  SYNDROME=", 
       TableForm[S]]]; Print["  "]; 
    NS = NullSpace[M, Modulus -> 2]; 
    ER = Sum[S[[i]]*2^(m - i), {i, m}]; 
    If[ER =!= 0, Print["POSITION OF ERROR=", ER]; Print[" "]; 
      W = w; W[[ER]] = Mod[W[[ER]] + 1, 2]; 
      Print["RECEIVED WORD=", w]; 
      Print["DECODED WORD= ", W]; Print[" "]; 
      Print["inf = ", LinearSolve[Transpose[NS], W, 
        Modulus -> 2]], W = w; 
      Print["THE RECEIVED WORD IS A CODE WORD"]; 
      Print["INFORM=", LinearSolve[Transpose[NS], W, 
        Modulus -> 2]]]]

ShowDecRectangularCode[n_Integer, m_Integer, w_List] := 
  Module[{s, A, d, OX, OY, er, i, j, DW, a, b, p}, 
   T = 0; If[TEST2[w] == 1, 
     Print["Check w on 0 or 1"]; Break[]]; 
    If[Length[w] != (n + 1)*(m + 1), 
     Print["Check the length of w"]; Break[]]; 
    Print["RECEIVED WORD=", w]; Print[" "]; 
    s = Partition[w, n + 1]; 
    Print["MATRIX FORM=", MatrixForm[s]]; Print[" "]; 
    OX = {}; OY = {}; d = 1; 
    Do[er = sum2[s[[i]]]; 
      If[er == 0, Print["ROW  ", i, " IS O.K."], 
       AppendTo[OX, i]], {i, m + 1}]; Print["  "]; 
    Do[er = sum2[Table[s[[i,j]], {i, m + 1}]]; 
      If[er == 0, Print["COLUMN ", j, " IS O.K."], 
       AppendTo[OY, j]], {j, n + 1}]; Print[" "]; 
    If[Length[OX] > 1 || Length[OY] > 1, 
     Print["2 OR MORE ERRORS"], 
     If[Length[OX] == 0, DW = s; d = 0, 
      p = s[[OX[[1]],OY[[1]]]]; 
       s[[OX[[1]],OY[[1]]]] = Mod[p + 1, 2]; DW = s; 
       Print["ERROR a[", OX[[1]], ",", OY[[1]], "]"]; d = 0]
      ]; Print["  "]; If[d == 1, Break[], 
     a = Drop[DW, {m + 1}]; b = Transpose[a]; 
      a = Drop[b, {n + 1}]; b = Transpose[a]; 
      a = Flatten[b]; Print["INFORM=", a]]]

ShowErrorTrappingDecoderBCHCode:=
  Module[{P,gr1,gr2,gr3,gr4,gr5,gr6,gr7,gr8,gr9,gr10,gr11,gr122,gr133,gr144,
      gr155,grtext,j},Do[P=Table["1",{8}];gr3=Table[Graphics[{}],{3}];
      gr7=Table[Graphics[{}],{2}];gr8=gr7;
      gr1=
        Graphics[{AbsoluteThickness[1],Line[{{-2,7},{63,7},{63,3}}],
            Line[{{0,
                  7},{0,-19},{19,-19},{19,-22},{59,-22},{59,-18},{19,-18},{19,
-20}}]}];
      gr2=
        Graphics[{AbsoluteThickness[1],Circle[{63,2},1],Text["+",{63,2}],
            Line[{{62,2},{1,
                  2},{1,-2},{1,-2},{3,-2},{3,-1},{5,-1},{5,-3},{3,-3},{3,-1}}]
,Line[{{32,2},{32,-1}}],Line[{{48,2},{48,-1}}],Line[{{56,2},{56,-1}}],
            If[j<16,Line[{{4,-3},{4,-5},{3,-6}}],Line[{{4,-3},{4,-8}}]]}];
      Do[gr3[[i]]=
          Graphics[{AbsoluteThickness[1],
              Line[{{5+8*(i-1),-2},{11+8*(i-1),-2},{11+8*(i-1),-1},{13+8*(i-1)
,-1},{13+8*(i-1),-3},{11+8*(i-1),-3},{11+8*(i-1),-1}}],
              If[j<16,Line[{{12+8*(i-1),-3},{12+8*(i-1),-5},{11+8*(i-1),-6}}],
                Line[{{12+8*(i-1),-3},{12+8*(i-1),-8}}]]}],{i,3}];
      gr4=Graphics[{AbsoluteThickness[1],Line[{{29,-2},{31,-2}}]}];
      gr5=
        Graphics[{AbsoluteThickness[1],Circle[{32,-2},1],Text["+",{32,-2}],
            Line[{{33,-2},{35,-2},{35,-1},{37,-1},{37,-3},{35,-3},{35,-1}}],
            If[j<16,Line[{{36,-3},{36,-5},{35,-6}}],
              Line[{{36,-3},{36,-8}}]]}];
      gr6=
        Graphics[{AbsoluteThickness[1],
            Line[{{37,-2},{43,-2},{43,-1},{45,-1},{45,-3},{43,-3},{43,-1}}],
            If[j<16,Line[{{44,-3},{44,-5},{43,-6}}],
              Line[{{44,-3},{44,-8}}]]}];
      Do[gr7[[i]]=
          Graphics[{AbsoluteThickness[1],
              Line[{{45+8*(i-1),-2},{47+8*(i-1),-2}}]}],{i,2}];
      Do[gr8[[i]]=
          Graphics[{AbsoluteThickness[1],Circle[{48+8*(i-1),-2},1],
              Text["+",{48+8*(i-1),-2}],
              Line[{{49+8*(i-1),-2},{51+8*(i-1),-2},{51+8*(i-1),-1},{53+8*(i-
1),-1},{53+8*(i-1),-3},{51+8*(i-1),-3},{51+8*(i-1),-1}}],
              If[i\[Equal]1,
                If[j<16,Line[{{52,-3},{52,-5},{51,-6}}],
                  Line[{{52,-3},{52,-8}}]],
                Line[{{60,-3},{60,-15},{62,-15}}]]}],{i,2}];
      gr9=Graphics[{AbsoluteThickness[1],Line[{{61,-2},{63,-2},{63,1}}]}];
      gr10=
        Graphics[{AbsoluteThickness[1],
            Line[{{1,-8},{55,-8},{55,-14},{1,-14},{1,-8}}],
            Line[{{28,-14},{28,-17},{62,-17}}]}];
      gr11=
        Graphics[{AbsoluteThickness[1],
            Line[{{62,-18},{62,-14},{65,-16},{62,-18}}],
            Line[{{65,-16},{67,-16},{67,2},{64,2}}],
            Line[{{67,-16},{67,-19}}]}];
      gr122=Graphics[{AbsoluteThickness[1],Line[{{59,-20},{66,-20}}]}];
      gr133=
        Graphics[{AbsoluteThickness[1],Circle[{67,-20},1],Text["+",{67,-20}],
            Line[{{68,-20},{70,-20}}],
            If[j<15||j>30,Line[{{70,-20},{70,-22},{71,-23}}],
              Line[{{70,-20},{70,-25},{10,-25},{10,-21},{19,-21}}]]}];
      gr144=
        Graphics[{Table[Line[{{4+8*(i-1),-6},{4+8*(i-1),-8}}],{i,7}],
            Line[{{70,-23},{70,-25},{10,-25},{10,-21},{19,-21}}],
            Line[{{73,-20},{75,-20}}]}];
      gr155=
        Graphics[{AbsoluteThickness[1],
            If[j>30,Line[{{70,-20},{75,-20}}],
              Line[{{70,-20},{72,-20},{73,-19}}]]}];
      grtext=
        Graphics[{Table[Text[P[[i]],{4+8*(i-1),-2}],{i,8}],Text["&",{63,-16}],
            Text[">",{30,7}],Text["<",{60,2}],Text["<",{10,2}],
            Text["v",{0,-5}],Text["<",{30,-25}],Text[">",{63,-20}],
            Text["^",{67,-5}],Text[">",{40,-17}],Text["out",{75,-20},{-1,0}],
            Text["v",{67,-17}],Text["in",{-3,7},{1,0}],Text["2",{30,-11}],
            Text["3",{35,-20}],
            If[j>15,Table[Text["v",{4+8*(i-1),-5}],{i,8}],Text[" ",{1,1}]]}];
      If[j\[Equal]1,Print[" "];Print["1 - the shift regesters"];Print[" "];
        Print["2 - the circuit of functional elements"];
        Print["    for the checking the syndrome on (000..0),"];
        Print["    or (000...010...0) only single 1 "];Print[" "];
        Print["3 - accumulative buffer"]];Print[" "];
      If[j\[Equal]1,Print["           1 The first 15 shifts."],
        If[j\[Equal]16,Print["           2 The next 15 shifts."],
          Print["           3 The last 15 shifts."]]];
      Show[gr1,gr2,gr3,gr4,gr5,gr6,gr7,gr8,gr9,gr10,gr11,gr122,gr133,gr144,
        gr155,grtext,PlotRange\[Rule]All],{j,1,31,15}]]

ShowHammingCode[m_Integer, inf_List] := 
  Module[{M, i, NS, W}, If[TEST2[inf] == 1, 
     Print["Check of the inf on 0 or 1"]; Break[]]; 
    If[Length[inf] != 2^m - 1 - m, 
     Print["Length of inf = ", Length[inf], "  =!=  ", 
       2^m - 1 - m]; Break[]]; 
    If[m > 4, Print["The parameter m is large"]; 
      Print["The parity check matrix will not be show on
the screen"]]; M = HammingMatrix[m]; 
    If[m <= 4, Print["PARITY"]; Print["CHECK"]; 
      Print["MATRIX=", MatrixForm[M]]]; Print["  "]; 
    NS = NullSpace[M, Modulus -> 2]; 
    If[m <= 4, Print["GENERATOR"]; 
      Print["MATRIX=", MatrixForm[NS]]; Print[" "]]; 
Print["inf = ", inf]; 
    Print["v = inf * G "]; W = Mod[inf . NS, 2]; 
    Print["HAMMING CODE VECTOR v =", W]]

ShowMeggittDecoder[g_,n_Integer,x_Symbol]:=
  Module[{cg,g1,k,i,j,l,gr1,gr2,gr3,gr4,gr6,gr7,gr8,gr9},
    g1=PolynomialMod[g,2];
    l=Length[CoefficientList[PolynomialMod[x^n-1,g1,Modulus\[Rule]2],{x}]];
    If[l\[NotEqual]0,
      Print["The polynomial g =",g1," isn`t the  generator polynomial,as"];
      Print["The generator polynomial of cyclic code of  length ",n,
        " must divides ",x^n-1];Break[]];Print[" "];
    cg=CoefficientList[g1,{x}];k=n-Length[cg]+1;
    gr1=Table[Graphics[PostScript[""]],{n-k}];
    gr2=Table[Graphics[PostScript[""]],{n-k}];
    gr3=Table[Graphics[PostScript[""]],{n-k}];
    gr1[[1]]=
      Graphics[{AbsoluteThickness[1],Circle[{-1,0},1],Text["+",{-1,0}],
          Text["^",{-2+8*(n-k),1.5}],
          Line[{{-2+8*(n-k),0},{-2+8*(n-k),2},{-1,2},{-1,1}}],
          Line[{{0,0},{2,0},{2,1},{4,1},{4,-1},{2,-1},{2,0}}],Text["<",{0,2}],
          Text[">",{1,0}],Text["1",{3,0}]}];
    Do[If[cg[[i]]\[Equal]1,
        gr2[[i]]=
          Graphics[{AbsoluteThickness[1],Line[{{-9+8*i,2},{-9+8*i,1}}],
              Text["v",{-9+8*i,1.5}]}]],{i,2,n-k}];
    Do[If[cg[[i]]\[Equal]1,
        gr1[[i]]=
          Graphics[{AbsoluteThickness[1],Text["+",{-9+8*i,0}],
              Circle[{-9+8*i,0},1],Text["",{-7+8*i,0}],Text["1",{-5+8*i,0}],
              Line[{{-8+8*i,0},{-6+8*i,0},{-6+8*i,1},{-4+8*i,
                    1},{-4+8*i,-1},{-6+8*i,-1},{-6+8*i,0}}]}],
        gr1[[i]]=
          Graphics[{AbsoluteThickness[1],Text["",{-7+8*i,0}],
              Line[{{-9+8*i,0},{-9+8*i,-3}}],Text["v",{-9+8*i,-2}],
              Text["1",{-5+8*i,0}],
              Line[{{-12+8*i,0},{-6+8*i,0},{-6+8*i,1},{-4+8*i,
                    1},{-4+8*i,-1},{-6+8*i,-1},{-6+8*i,0}}]}]],{i,2,n-k}];
    Do[If[cg[[i+1]]\[Equal]1,
        gr3[[i]]=
          Graphics[{AbsoluteThickness[1],Line[{{-4+8*i,0},{-2+8*i,0}}],
              Line[{{-3+8*i,0},{-3+8*i,-3}}],Text["v",{-3+8*i,-2}],
              Text["",{-2.5+8*i,0}]}]],{i,n-k}];
    gr4=Graphics[{AbsoluteThickness[1],Line[{{-10,-10},{-5,-10}}],
          Line[{{-2,0},{-5,
                0},{-5,-10},{-2,-10},{-2,-11},{-2+8*(n-k),-11},{-2+8*(n-
                        k),-9},{-2,-9},{-2,-10}}],
          Line[{{-2+8*(n-k),-10},{2+8*(n-k),-10}}]}];
    gr6=Graphics[{AbsoluteThickness[1],Circle[{3+8*(n-k),-10},1],
          Text["+",{3+8*(n-k),-10}],
          Line[{{4+8*(n-k),-10},{6+8*(n-k),-10}}]}];
    gr7=Graphics[{AbsoluteThickness[1],
          Line[{{2,-3},{-2+8*(n-k),-3},{-2+8*(n-k),-7},{2,-7},{2,-3}}],
          Line[{{-2+8*(n-k),-5},{3+8*(n-k),-5},{3+8*(n-k),-9}}]}];
    gr8=Graphics[{Text["^",{-5,-5}],Text[">",{-3,-10}],
          Text[">",{8*(n-k),-10}],Text[">",{5+8*(n-k),-10}],
          Text["v",{3+8*(n-k),-7}]}];
    If[n-k>2,Print["1 - the shift register"];Print[" "];
      gr9=
        Graphics[{Text["Accumulative buffer",{-2+4*(n-k),-10}],
            Text["The circuit of",{4*(n-k),-3.5}],
            Text["functional ele-",{4*(n-k),-4.5}],
            Text["ments for the",{4*(n-k),-5.5}],
            Text["search of errors",{4*(n-k),-6.5}]}],
      gr9=Graphics[{Text["3",{-2+4*(n-k),-10}],Text["2",{4*(n-k),-5}]}];
      Print["1 - the shift register"];
      Print[" "];Print["2 - the circuit of functional elements"];
      Print["    for the search of errors"];Print[" "];
      Print["3 - the accumulative buffer"];Print[" "]];
    Show[gr1,gr2,gr3,gr4,gr6,gr7,gr8,gr9,PlotRange\[Rule]All]]

ShowNonsystematicEncoderCyclicCode[g_,n_Integer,x_Symbol]:=
  Module[{g1,cg,i,j,k,ph,gr,A,B},g1=PolynomialMod[g,2];
    l=Length[CoefficientList[PolynomialMod[x^n-1,g1,Modulus\[Rule]2],{x}]];
    If[l\[NotEqual]0,
      Print["The polynomial g= ",g1," isn`t a generator polynomial, as"];
      Print["the generator polynomial of cyclic code of length ",n,
        " must divides ",x^n-1];Break[]];
    cg=Reverse[CoefficientList[g1,{x}]];k=Length[cg];
    Print["GENERATOR POLYNOMIAL = ",g1];gr=Table["",{k}];
    ph=Table["",{k-1}];
    gr[[1]]=Graphics[{AbsoluteThickness[1],Text[">",{-9,0}],
          Line[{{-11,0},{-3,0},{-3,4},{3,4}}],Text["^",{-3,2}],
          Line[{{-3,0},{1,0}}],Text[">",{0,0}],
          Line[{{1,1},{3,1},{3,-1},{1,-1},{1,1}}]}];
    Do[gr[[i]]=
        If[cg[[i]]\[Equal]0,
          Graphics[{AbsoluteThickness[1],Line[{{4*i-5,0},{4*i-3,0}}],
              Text[">",{4*i-3.5,0}],
              Line[{{4*i-3,1},{4*i-1,1},{4*i-1,-1},{4*i-3,-1},{4*i-3,1}}]}],
          Graphics[{AbsoluteThickness[1],Line[{{4*i-5,0},{4*i-3,0}}],
              Text[">",{4*i-3.5,0}],
              Line[{{4*i-3,1},{4*i-1,1},{4*i-1,-1},{4*i-3,-1},{4*i-3,1}}],
              Text["^",{4*i-4,2}],Line[{{4*i-4,0},{4*i-4,3}}]}]],{i,2,k-1}];
    gr[[k]]=
      Graphics[{AbsoluteThickness[1],Text[">",{4*k-2,0}],Text["^",{4*k-1,2}],
          Line[{{4*(k-1)-1,0},{4*(k-1)+3,0},{4*(k-1)+3,3}}]}];
    Do[ph[[i]]=
        If[cg[[i+1]]\[Equal]0,
          Graphics[{AbsoluteThickness[1],Line[{{4*i-1,4},{4*i+3,4}}]}],
          Graphics[{AbsoluteThickness[1],Text[">",{4*i-1.5,4}],
              Circle[{4*i,4},1],Text["+",{4*i,4}],
              Line[{{4*i+1,4},{4*i+3,4}}]}]],{i,k-2}];
    ph[[k-1]]=
      Graphics[{AbsoluteThickness[1],Text[">",{4*k-2.5,4}],
          Text[">",{4*k+7,4}],Circle[{4*(k-1)+3,4},1],Text["+",{4*(k-1)+3,4}],
          Line[{{4*k,4},{4*k+10,4}}],Line[{{4*k-2,4},{4*k-5.1,4}}]}];
    A=Graphics[{Text["inf",{-4,-0.5},{1,0}]}];
    B=Graphics[{Text["code word",{4*k+3,4.5},{-1,0}]}];
    Show[ph,gr,A,B,PlotRange\[Rule]All]]


ShowReconstructInformRectangularCode[n_Integer, m_Integer, 
   w_List,X_Symbol] := 
  Module[{i, j, DW, a, b}, 
   If[Length[w] != (n + 1)*(m + 1), 
     Print["CHECK THE LENGTH w =          ", 
       Length[w], " != ", (n + 1)*(m + 1)]; 
      Break[]]; If[TEST2[w] == 1, 
     Print["CHECKw ON 0 OR 1"]; Break[]]; 
    Print["THE CORRECTING WORD=", w]; Print[" "]; 
    DW = Table[0, {m + 1}, {n + 1}]; 
    Do[Do[DW[[i,j]] = 
       w[[Mod[(m + 1)*(j - i) + i - 1, 
          (m + 1)*(n + 1)] + 1]], {i, m + 1}], 
     {j, n + 1}]; Print["MATRIX FORM=", 
     MatrixForm[DW]]; Print[" "]; 
    a = Drop[DW, {m + 1}]; b = Transpose[a]; 
    a = Drop[b, {n + 1}]; b = Transpose[a]; 
    a = Flatten[b]; Print["INFORM=", a]]

ShowRectangularCode[n_Integer, m_Integer, inf_List,X_Symbol] := 
  Module[{a, i, j, w, b, c}, 
   If[TEST2[inf] == 1, Print["CHECK inf on 0 or 1 "]; 
      Break[]]; If[Length[inf] != n*m, 
     Print["CHECK THE LENGTH inf"]; Break[]]; a = inf; 
    Print["INFORM=", a]; Print["  "]; 
    Do[a = Insert[a, X, i*(n + 1)], {i, m}]; 
    Do[AppendTo[a, X], {i, n + 1}]; a = Partition[a, n + 1]; 
    Print["MATRIX FORM=", MatrixForm[a]]; Print["  "]; 
    Do[a[[i,n + 1]] = sum2[Drop[a[[i]], {n + 1}]], {i, m}]; 
    Do[a[[m + 1,j]] = sum2[Table[a[[i,j]], {i, m}]], 
     {j, n + 1}]; Print["MATRIX FORM=", MatrixForm[a]]; 
    Print[" "]; w = Flatten[a]; Print["CODE WORD=", w]]

ShowStandardArrays[M_List, w_List] := 
  Module[{st, p, NS, n, i, j, m, A, B, inf, k, K, cd, lcs, 
    W, r, ST, c}, If[TEST1[M] == 1, 
     Print["Check the matrix M "]; Break[]]; 
    If[TEST2[M] == 1, Print["Check the matrix M on 0 or 1"]; 
      Break[]]; If[TEST2[w] == 1, 
     Print["Check vactor w on 0 or 1"]; Break[]]; 
    If[Length[w] != Length[M[[1]]], 
     Print["Check the length of w"]; Break[]]; Print[" "]; 
    NS = NullSpace[M, Modulus -> 2]; n = Length[NS]; 
    m = Length[NS[[1]]]; 
    A = Array[X, {2^(m - n), 2^n}]; 
    B = Table[VEC[i, m], {i, 0, 2^m - 1}]; 
    inf = Table[VEC[i, n], {i, 0, 2^n - 1}]; 
    Print["PARITY CHECK MATRIX=", MatrixForm[M]]; 
    A[[1]] = Table[Mod[Sum[inf[[i,j]]*NS[[j]], {j, n}], 2], 
      {i, 2^n}]; Do[Do[k = Position[B, A[[i,j]]]; 
        K = k[[1,1]]; B = Drop[B, {K, K}], {j, 2^n}]; 
      cd = Table[Mod[B[[1]] + A[[1,j]], 2], {j, 2^n}]; 
      lcs = cd[[1]]; Do[If[Apply[Plus, cd[[i]]] < 
         Apply[Plus, lcs], lcs = cd[[i]], Continue[]], 
       {i, 2, 2^n}]; A[[i + 1]] = 
       Table[Mod[lcs + A[[1,j]], 2], {j, 2^n}], 
     {i, 2^(m - n) - 1}]; 
    st = Table[" ", {2^(m - n)}, {2^n}]; 
    Do[Do[Do[If[A[[i,j]][[p]] == 0, 
        st[[i,j]] = StringInsert[st[[i,j]], "0", -1], 
        st[[i,j]] = StringInsert[st[[i,j]], "1", -1]], 
       {p, m}], {i, 2^(m - n)}], {j, 2^n}]; W = " "; 
    Do[If[w[[i]] == 0, W = StringInsert[W, "0", -1], 
      W = StringInsert[W, "1", -1]], {i, Length[w]}]; 
    r = StringPosition[StringJoin[st], W]; 
    c = (r[[1,1]] - 1)/StringLength[W] + 1; 
    If[Mod[c, 2^n] == 0, j = 2^n, j = Mod[c, 2^n]]; 
    i = (c - j)/2^n + 1; If[i == 0, i = 1]; ST = st; 
    ST[[i,j]] = " "; Do[ST[[i,j]] = 
      StringInsert[ST[[i,j]], "X", -1], {k, m}]; 
    If[m > 6 || 2^n > 8, 
     Print["RECEIVED WORD =", W]; Print[" "]; 
      Print[" Code             Coset"]; Print[" Words"]; 
      Print[TableForm[st[[1]]], "        ", 
       TableForm[ST[[i]]], "<-- Coset leader"], 
     Print["Coset"]; Print["Leader"]; 
      Print[MatrixForm[st]]; Print[" "]; 
      Print["RECEIVED WORD =", W]; Print[" "]; 
      Print[MatrixForm[ST]]]; Print[" "]; 
    Print["DECODED WORD=", st[[1,j]]]]

ShowSystematicEncode[M_List, inf_List, x_Symbol] := 
  Module[{v, u, a, n1, l, d, i, S, w, A, B, H, p, k, g, j, vi}, 
    If[TEST1[M] == 1, Print["CHECK THE MATRIX M "]; Break[]];
    If[TEST2[M] == 1, Print["CHECK M ON 0 OR 1"]; Break[]];
    If[TEST2[inf] == 1, Print["CHECK inf ON 0 OR 1"]; Break[]];
    If[DimensionCode[M] != Length[inf], 
      Print["Check the length of information word"]; Break[]];
    r = 0;
    If[Length[M[[1]]] > 18, Print["The parametrs of matrix H are large,"];
      Print["the parity  check matrix will not be show on  the screen"]; 
      r = 1]; Print[" "];
    If[r == 0, Print["Parity Check Matrix =", M]; 
      v = Array[x, {Length[M[[1]]]}];
      Print["The code word =", v]; u = Table[0, {Length[M]}];
      l = 
        Table[v[[i]], {i, Length[M[[1]]] - Length[inf] + 1, Length[M[[1]]]}];
      a = M.v == u;
      S = Solve[a && Modulus == 2, l];
      S = Flatten[S];
      S = Rest[S];
      w = v /. S;
      vi = 
        Flatten[Transpose[
            Position[Table[w[[i]] == v[[i]], {i, 1, Length[w]}], True]]];
      If[vi == Table[i, {i, 1, Length[M[[1]]] - Length[S]}], n1 = 1, n1 = 0];
      If[n1 == 1, d = Table[v[[i]], {i, Length[M[[1]]] - Length[S]}];
        l = 
          Table[v[[i]], {i, Length[M[[1]]] - Length[S] + 1, Length[M[[1]]]}];
        Print["The informations symbols =", d];
        Print["The check symbols = ", l];
        Print[" "];
        Print["The check symbols compute using this equations"];
        Print[TableForm[S]]; Print[" "]; Print["The code word=", w]; 
        Print[" "];
        a = Table[d[[i]] -> inf[[i]], {i, Length[d]}];
        Print["The inf=      ", inf];
        Print["The code word=", Mod[w /. a, 2]] , Print[TableForm[S]];
        Print["This code is not  systematic, but, if you need", 
          "the equivalent systematic code K1 you must do a permutation", 
          "of numbers  1,2,...such that for each word v  we", 
          "have:(v1,v2,...,vn) in the code K <=>", 
          "(vg(1),vg(2),...,vg(n)) in the code K1"];
        g = Join[vi, Complement[Table[i, {i, 1, Length[M[[1]]]}], vi]];
        Print["g =", g];
        Print[" "];
        Print["If you want continue the calculation"];
        Print["you must try"];
        Print["M=",M1 = EquivalentCode[M, g],"; ShowSystematicEncode[M,inf,x]"]]]]

ShowSystematicEncoderCyclicCode[g_,n_Integer,x_Symbol]:=
  Module[{cg,k,gr1,gr2,gr3,gr4,gr5,gr6,i,j,l},g1=PolynomialMod[g,2];
    l=Length[CoefficientList[PolynomialMod[x^n-1,g1,Modulus\[Rule]2],{x}]];
    If[l\[NotEqual]0,
      Print["The polynomial g =",g1," isn`t the  generator polynomial,as"];
      Print["The generator polynomial of cyclic code of  length ",n,
        " must divides ",x^n-1];Break[]];Print[" "];
    cg=CoefficientList[g1,{x}];k=n-Length[cg]+1;
    gr1=Table["",{n-k}];gr2=Table[Graphics[PostScript[""]],{n-k}];
    gr3=Table[Graphics[PostScript[""]],{n-k}];
    Print["              1. The first ",k," shifts."];
    Print["        (the length of information word =",k,")"];
    gr1[[1]]=
      Graphics[{AbsoluteThickness[1],Text["+",{-9+8*(n-k+1),0}],
          Circle[{-9+8*(n-k+1),0},1],Line[{{-1+8*(n-k),1},{-1+8*(n-k),3}}],
          Text["^",{-1+8*(n-k),2}],
          Line[{{-6+8*(n-k),3},{-1,3},{-1,0},{2,0},{2,1},{4,
                1},{4,-1},{2,-1},{2,0}}],Text["<",{0,3}],Text[">",{1,0}]}];
    Do[If[cg[[i]]\[Equal]1,
        gr2[[i]]=
          Graphics[{AbsoluteThickness[1],Line[{{-9+8*i,3},{-9+8*i,1}}],
              Text["v",{-9+8*i,2}]}]],{i,n-k}];
    Do[If[cg[[i]]\[Equal]1,
        gr1[[i]]=
          Graphics[{AbsoluteThickness[1],Text["+",{-9+8*i,0}],
              Circle[{-9+8*i,0},1],Text[">",{-7+8*i,0}],
              Line[{{-8+8*i,0},{-6+8*i,0},{-6+8*i,1},{-4+8*i,
                    1},{-4+8*i,-1},{-6+8*i,-1},{-6+8*i,0}}]}],
        gr1[[i]]=
          Graphics[{AbsoluteThickness[1],Text[">",{-7+8*i,0}],
              Line[{{-9+8*i,0},{-6+8*i,0},{-6+8*i,1},{-4+8*i,
                    1},{-4+8*i,-1},{-6+8*i,-1},{-6+8*i,0}}]}]],{i,2,n-k}];
    Do[If[cg[[i+1]]\[Equal]1,
        gr3[[i]]=
          Graphics[{AbsoluteThickness[1],Line[{{-4+8*i,0},{-2+8*i,0}}],
              Text[">",{-2.5+8*i,0}]}],
        gr3[[i]]=
          Graphics[{AbsoluteThickness[1],Line[{{-4+8*i,0},{2+8*i,0}}]}],
        Text[">",{-2.5+8*i,0}]],{i,n-k}];
    gr4=Graphics[{AbsoluteThickness[1],
          Text["The information word->",{-1,-2.5},{-1,0}],Text[">",{0,-2}],
          Line[{{-1,-2},{2+8*(n-k),-2},{-1+8*(n-k),-2},{-1+8*(n-k),-1}}],
          Text["^",{-1+8*(n-k),-1.5}]}];
    gr5=Graphics[{AbsoluteThickness[1],Line[{{8*(n-k),0},{8*(n-k)+2,0}}],
          Text["The code",{8*(n-k)+3.7,-0.7},{-1,0}],
          Text["word",{8*(n-k)+6,-1.2},{-1,0}]}];
    Do[If[j\[Equal]1,
        gr6=Graphics[{AbsoluteThickness[1],
              Line[{{8*(n-k)+2,-2},{8*(n-k)+4,-1},{8*(n-k)+10,-1}}],
              Line[{{-1+8*(n-k),3},{-6+8*(n-k),3}}],
              Text[">",{8*(n-k)+5,-1}]}],
        gr6=Graphics[{AbsoluteThickness[1],Text[">",{8*(n-k)+5,-1}],
              Line[{{8*(n-k)+2,0},{8*(n-k)+4,-1},{8*(n-k)+10,-1}}],
              Line[{{-1+8*(n-k),3},{-5+8*(n-k),3.5}}]}]];Print[" "];
      Show[gr1,gr2,gr3,gr4,gr5,gr6,PlotRange\[Rule]All];Print[" "];
      If[j\[Equal]1,Print["              2. The last ",n-k," shifts"]],{j,
        2}]]

ShowTableMult[p_Integer] := 
  Module[{T, i, j}, T = Table[0, {p + 2}, {p + 2}]; 
    T[[1,1]] = "*"; Do[T[[1,i]] = i - 3; T[[i,1]] = i - 3, 
     {i, 3, p + 2}]; Do[T[[2,i]] = ""; T[[i,2]] = "", 
     {i, p + 2}]; Do[Do[T[[i,j]] = 
       Mod[(i - 3)*(j - 3), p], {i, 3, p + 2}], 
     {j, 3, p + 2}]; Print[" "]; Print[MatrixForm[T]]]

ShowTablePlus[p_Integer] := 
  Module[{T, i, j}, T = Table[0, {p + 2}, {p + 2}]; 
    T[[1,1]] = "+"; Do[T[[1,i]] = i - 3; T[[i,1]] = i - 3, 
     {i, 3, p + 2}]; Do[T[[2,i]] = ""; T[[i,2]] = "", 
     {i, p + 2}]; Do[Do[T[[i,j]] = Mod[i + j - 6, p], 
      {i, 3, p + 2}], {j, 3, p + 2}]; Print[" "]; 
    Print[MatrixForm[T]]]

SingletonBound[n_,k_,d_]:=Module[{},
If[IntegerQ[n]&&IntegerQ[k]&&IntegerQ[d],
If[n-k-d+1>=0,True,False],
 If[IntegerQ[n]&&IntegerQ[k],d <= n-k+1,
  If[IntegerQ[n]&&IntegerQ[d],k <=n-d+1,
   If[IntegerQ[k]&&IntegerQ[d],n >= k+d-1,
    If[IntegerQ[d],n >= k+d-1,
     If[IntegerQ[k],n >=d+ k-1,
      If[IntegerQ[n],k + d <= n+1,
n - k >= d - 1]]]]]]]]

Syndrome[H_List,w_List]:= Mod[H.w,2]

SystematicEncode[M_List, inf_List,x_Symbol] := 
  Module[{v, u, a, n1, l, d, i, S, w, A, B, H, p, k, g, j, vi,cw}, 
    If[TEST1[M] == 1, Print["CHECK THE MATRIX M "]; Break[]];
    If[TEST2[M] == 1, Print["CHECK M ON 0 OR 1"]; Break[]];
    If[TEST2[inf] == 1, Print["CHECK inf ON 0 OR 1"]; Break[]];
    If[DimensionCode[M] != Length[inf], 
      Print["Check the length of information word"]; Break[]];
      v = Array[x, {Length[M[[1]]]}];u = Table[0, {Length[M]}];
      l = 
        Table[v[[i]], {i, Length[M[[1]]] - Length[inf] + 1, Length[M[[1]]]}];
      a = M.v == u;
      S = Solve[a && Modulus == 2, l];
      S = Flatten[S];
      S = Rest[S];
      w = v /. S;
      vi = 
        Flatten[Transpose[
            Position[Table[w[[i]] == v[[i]], {i, 1, Length[w]}], True]]];
      If[vi == Table[i, {i, 1, Length[M[[1]]] - Length[S]}], n1 = 1, n1 = 0];
      If[n1 == 1, d = Table[v[[i]], {i, Length[M[[1]]] - Length[S]}];
        l = Table[v[[i]], {i, Length[M[[1]]] - Length[S] + 1, Length[M[[1]]]}];
        a = Table[d[[i]] -> inf[[i]], {i, Length[d]}];
        cw=Mod[w /. a, 2] , 
        Print["This code is not  systematic, but, if you need", 
          "the equivalent systematic code K1 you must do a permutation", 
          "of numbers  1,2,...such that for each word v  we", 
          "have:(v1,v2,...,vn) in the code K <=>", 
          "(vg(1),vg(2),...,vg(n)) in the code K1"];
        g = Join[vi, Complement[Table[i, {i, 1, Length[M[[1]]]}], vi]];
        Print["g =", g];
        Print[" "];
        Print["If you want continue the calculation"];
        Print["you must try"];
        Print["M=",M1 = EquivalentCode[M, g],"; SystematicEncode[M,inf,x]"];If[n1==1,cw]]]


SystematicEncodeCyclicCode[g_,infword_List,n_Integer,x_Symbol]:=
  Module[{cg,k,gr1,gr2,gr3,gr4,gr5,gr6,gr7,gr8,i,j,l},inf=Reverse[infword];
    g1=PolynomialMod[g,2];
    l=Length[CoefficientList[PolynomialMod[x^n-1,g1,Modulus\[Rule]2],{x}]];
    If[l\[NotEqual]0,
      Print["The polynomial g =",g1," isn`t the  generator polynomial,as"];
      Print["The generator polynomial of cyclic code of  length ",n,
        " must divides ",x^n-1];Break[]];Print[" "];
    cg=CoefficientList[g1,{x}];k=n-Length[cg]+1;
    gr1=Table["",{n-k}];gr2=Table[Graphics[PostScript[""]],{n-k}];
    gr3=Table[Graphics[PostScript[""]],{n-k}];gr8=Graphics[PostScript[""]];
    If[TEST2[inf]\[Equal]1,Print["Check inf on 0 or 1."];Break[]];
    If[Length[inf]\[NotEqual]k,
      Print["The length inf =",Length[inf],",but != ",k,"."];Break[]];
    pl=Sum[inf[[i]]*x^(n-k+i-1),{i,k}];
    Print["THE INFORMATION POLYNOMIAL = ",pl];L=inf;
    Do[PrependTo[L,0],{n-k}];RL=Reverse[L];L1="";
    Do[If[inf[[i]]\[Equal]0,L1=StringInsert[L1,"0",-1],
        L1=StringInsert[L1,"1",-1]],{i,k}];
    Do[L1=StringInsert[L1," ",1],{n-k}];
    cpl=pl+PolynomialMod[pl,g,Modulus\[Rule]2];cw="";j1=1;
    cdw=CoefficientList[cpl,{x}];r=Length[cdw];
    Do[AppendTo[cdw,0],{n-r}];
    Do[If[cdw[[i]]\[Equal]0,cw=StringInsert[cw,"0",-1],
        cw=StringInsert[cw,"1",-1]],{i,n}];P=Table["A",{n-k}];
    NGr1=Table[0,{n},{n-k}];NGr4=RL;NGr6=Reverse[cdw];
    NGr7=Table[0,{n}];W=P;
    Do[NGr1[[j,1]]=If[IntegerQ[P[[n-k]]],Mod[RL[[j]]+P[[n-k]],2],RL[[j]]];
      If[j\[LessEqual]k,W[[1]]=NGr1[[j,1]];
        gr1[[1]]=
          Graphics[{AbsoluteThickness[1+NGr1[[j,1]]],
              Text["+",{-9+8*(n-k+1),0}],Circle[{-9+8*(n-k+1),0},1],
              Line[{{-1+8*(n-k),1},{-1+8*(n-k),3}}],Text["^",{-1+8*(n-k),2}],
              Line[{{-6+8*(n-k),3},{-1,3},{-1,0},{2,0},{2,1},{4,
                    1},{4,-1},{2,-1},{2,0}}],Text[W[[1]],{3,0}],
              Text["<",{0,3}],Text[">",{1,0}]}],W[[1]]=" ";
        gr1[[1]]=
          Graphics[{AbsoluteThickness[1+NGr1[[j,1]]],
              Text["+",{-9+8*(n-k+1),0}],Circle[{-9+8*(n-k+1),0},1],
              Line[{{-1+8*(n-k),1},{-1+8*(n-k),3}}],
              Text["^",{-1+8*(n-k),2}]}];
        gr8=
          Graphics[{AbsoluteThickness[1],Text[W[[1]],{3,0}],
              Line[{{-6+8*(n-k),3},{-1,3},{-1,0},{2,0},{2,1},{4,
                    1},{4,-1},{2,-1},{2,0}}]}]];
      Do[If[cg[[i]]\[Equal]1,
          gr2[[i]]=
            Graphics[{AbsoluteThickness[1+If[j\[LessEqual]k,NGr1[[j,1]],0]],
                Line[{{-9+8*i,3},{-9+8*i,1}}],Text["v",{-9+8*i,2}]}]],{i,
          n-k}];Do[
        If[cg[[i]]\[Equal]1,
          W[[i]]=If[j\[LessEqual]k,
              Mod[NGr1[[j,1]]+If[IntegerQ[P[[i-1]]],P[[i-1]],0],2],
              If[IntegerQ[P[[i-1]]],P[[i-1]],""]];
          
          gr1[[i]]=
            Graphics[{AbsoluteThickness[
                  1+If[j\[LessEqual]k,W[[i]],
                      If[IntegerQ[P[[i-1]]],P[[i-1]],0]]],
                Text[W[[i]],{-5+8*i,0}],Text["+",{-9+8*i,0}],
                Circle[{-9+8*i,0},1],Text[">",{-7+8*i,0}],
                Line[{{-8+8*i,0},{-6+8*i,0},{-6+8*i,1},{-4+8*i,
                      1},{-4+8*i,-1},{-6+8*i,-1},{-6+8*i,0}}]}],
          W[[i]]=If[IntegerQ[P[[i-1]]],P[[i-1]],""];
          
          gr1[[i]]=
            Graphics[{AbsoluteThickness[1+If[IntegerQ[P[[i-1]]],P[[i-1]],0]],
                Text[">",{-7+8*i,0}],Text[W[[i]],{-5+8*i,0}],
                Line[{{-12+8*i,0},{-6+8*i,0},{-6+8*i,1},{-4+8*i,
                      1},{-4+8*i,-1},{-6+8*i,-1},{-6+8*i,0}}]}]],{i,2,n-k}];
      Do[If[cg[[i+1]]\[Equal]1,
          gr3[[i]]=
            Graphics[{AbsoluteThickness[1+If[IntegerQ[P[[i]]],P[[i]],0]],
                Line[{{-4+8*i,0},{-2+8*i,0}}],Text[">",{-2.5+8*i,0}]}]],{i,
          n-k}];gr4=
        Graphics[{AbsoluteThickness[1+RL[[j]]],
            Text[StringTake[StringReverse[L1],{j}],{8*(n-k)-4,-2}],
            Text[StringTake[L1,n-j],{-1,-2.5},{-1,0}],Text[">",{0,-2}],
            Line[{{-1,-2},{8*(n-k)-5,-2}}],
            Line[{{8*(n-k)-3,-2},{2+8*(n-k),-2},{-1+8*(n-k),-2},{-1+8*(n-
                          k),-1}}],Text["^",{-1+8*(n-k),-1.5}]}];
      gr5=
        Graphics[{AbsoluteThickness[1+NGr1[[j,1]]],
            Line[{{8*(n-k),0},{8*(n-k)+2,0}}],
            Text[StringTake[cw,1-j],{8*(n-k)+7,-0.7},{-1,0}]}];
      If[j>k,j1=0];
      If[j1\[Equal]1,
        gr6=Graphics[{AbsoluteThickness[1+NGr6[[j]]],
              Line[{{8*(n-k)+2,-2},{8*(n-k)+4,-1},{8*(n-k)+5,-1}}],
              Line[{{8*(n-k)+7,-1},{8*(n-k)+12,-1}}]}];
        gr7=
          Graphics[{AbsoluteThickness[1+NGr1[[j,1]]],
              Line[{{-1+8*(n-k),3},{-6+8*(n-k),3}}],
              Text[StringTake[cw,{-j}],{8*(n-k)+6,-1}]}],
        gr6=Graphics[{AbsoluteThickness[1+NGr6[[j]]],
              Text[StringTake[cw,{-j}],{8*(n-k)+6,-1}],
              Line[{{8*(n-k)+2,0},{8*(n-k)+4,-1},{8*(n-k)+5,-1}}],
              Line[{{8*(n-k)+7,-1},{8*(n-k)+12,-1}}]}];
        gr7=
          Graphics[{AbsoluteThickness[1+NGr1[[j,1]]],
              Line[{{-1+8*(n-k),3},{-5+8*(n-k),3.5}}]}]];Print[" "];
      Show[gr1,gr2,gr3,gr4,gr5,gr6,gr7,gr8];P=W,{j,n}];
    Print["THE CODE POLYNOMIAL = ",cpl];Print["THE CODE WORD = ",StringReverse[cw]]]


VEC[n_, m_] := Module[{a, l, i}, If[n>2^m,Print[n," > ",2^m];Break[]];
   i = n; If[i == 0, a = {0}, a = {}; 
      While[i > 0, If[IntegerQ[i/2], PrependTo[a, 0]; i /= 2, 
        PrependTo[a, 1]; i = (i - 1)/2]]]; l = Length[a]; 
    Do[PrependTo[a, 0], {i, m - l}]; a]

WeightPolynomial[M_List,x_Symbol] := 
  Module[{a, i, A, S, m, n, k, NS, b}, 
   If[TEST1[M] == 1, Print["Check the matrix M"]; Break[]]; 
    If[TEST2[M] == 1, Print["Check the matrix on 0 or 1"]; 
      Break[]];
    If[DimensionCode[M] <=DimensionCode[DualCode[M] ], 
      NS = NullSpace[M, Modulus -> 2]; 
    n = Length[NS[[1]]]; m = Length[NS]; 
    A = Table[0, {n + 1}]; 
    Do[a = VEC[i, m]; k = 
       Apply[Plus,Mod[a.NS, 2]]; 
      A[[k + 1]] += 1, {i, 0, 2^m - 1}]; 
    S = Sum[A[[i + 1]]*x^i, {i, 0, n}]; S,
      NS = NullSpace[DualCode[M] , Modulus -> 2]; 
    n = Length[NS[[1]]]; m = Length[NS]; 
    A = Table[0, {n + 1}]; 
    Do[a = VEC[i, m]; k = 
       Apply[Plus,Mod[a.NS, 2]]; 
      A[[k + 1]] += 1, {i, 0, 2^m - 1}]; 
    S = Sum[A[[i + 1]]*x^i, {i, 0, n}]; S=McWilliamsIdentity[S,n,x] ];S]
  

(*-------------------------------------*)

SYNDROMLIST[g_, n_Integer, k_Integer, d_Integer,x_Symbol] := 
  Module[{j, c, L, l, i}, If[d > 8, 
     Print["The parameters of code are large,"]; 
      Print["try decode the another method."]; Clear[Distance, Dsure, Stop]; 
      Break[]]; j = Floor[(d - 1)/2.] - 1; dis = d; Clear[Distance, Dsure]; 
    t = j + 1; If[j == -1, Print["This code cannot correct"]; 
      Print["even a single error"]; Clear[Distance, Stop, Dsure]; Break[]]; 
    L = Sum[Binomial[n - 1, i], {i, 0, j}]; c = Table[0, {L}]; 
    Do[c[[l]] = 
       CoefficientList[PolynomialMod[x^(n - 1) + 
          If[l <= n, If[l == 1, 0, 1]*x^(l - 2), 
           Sum[First[Take[Flatten[Take[Permutations[Table[If[i <= 
                     n - 3, 0, 1], {i, n - 1}]], {l - n, l - n}]], 
               {i, i}]]*x^(i - 1), {i, n - 1}]], g, Modulus -> 2], {x}]; 
      Do[AppendTo[c[[l]], 0], {n - k - Length[c[[l]]]}]; 
      If[Mod[l, 50] == 0, Print[l, " Please wait I find the Syndrome list"]], 
     {l, L}]; If[Length[c] != Length[Union[c]], 
     Print["This code has another a code distance"]; 
      Print["if you want, you can change a code distance and try again,"]; 
      Print["Distance =  ; MeggittDecode[g,w,n,x]"]; Break[]]; c]

TEST1[M_] := 
  Module[{i, j, T}, T = 0; 
    Do[Do[If[Length[M[[j]]] != Length[M[[i]]], 
       T = 1; Break[]], {i, j + 1, Length[M]}], 
     {j, Length[M]}]; T]


TEST2[M_] := 
  Module[{H, T, i}, T = 0; H = Flatten[M]; 
    Do[If[H[[i]] != 0 && H[[i]] != 1, T = 1; Break[]]; 
      If[H[[i]] == Null, T = 1; Break[]]; 
      If[IntegerQ[H[[i]]], Continue[], T = 1; Break[]], 
     {i, Length[H]}]; T]

sum2[a_] := Mod[Apply[Plus, a], 2]

POLPRINTLENGTH[g_,x_Symbol] := 
  Module[{l}, If[IntegerQ[g], 1, 
    l = CoefficientList[g, {x}]; 
     If[l[[1]] == 1, If[l[[2]] == 0, 
       5*(g /. x -> 1) - 4, 5*(g /. x -> 1) - 5], 
      If[l[[2]] == 0, 5*(g /. x -> 1) - 3, 
       5*(g /. x -> 1) - 4]]]]

PP[k_, n_] := If[Mod[k, n] == 0, n, Mod[k, n]]

End[ ]
Print[" ready.."]
EndPackage[ ]