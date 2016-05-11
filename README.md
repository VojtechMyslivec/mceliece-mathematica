# Asymetrický šifrovací algoritmus McEliece

**Diplomová práce**

 - České vysoké učení technické v Praze
 - Fakulta informačních technologií
 - Letní semestr 2016


## Abstrakt

Česká verze níže

### English

In this work, we deal with a code-based public-key cryptosystem *McEliece* which
is one of the candidates for post-quantum cryptography. We provide a definition
of the cryptosystem, its variant for digital signature scheme, and we focus on
the practical aspects of this cryptosystem and its cryptanalysis. We evaluate
the time complexity of the algorithms using an illustrative implementation in
*Wolfram Mathematica*.

### Česky

V této práci se zabýváme asymetrickým kryptosystémem *McEliece*, který je
založený na samoopravných lineárních kódech a je jedním z kandidátů pro
asymetrickou postkvantovou kryptografii. V práci uvádíme základní definici
tohoto kryptosystému, variantu pro digitální podpis a věnujeme se též
existujícím kryptoanalýzám a praktickým aspektům tohoto systému. V rámci práce
vznikla ukázková implementace v softwaru *Wolfram Mathematica*, na které bylo
provedeno měření časových závislostí algoritmů.

## Zadání

Prostudujte asymetrický šifrovací algoritmus McEliece založený na binárních
Goppa kódech. Proveďte rešerši existujících kryptoanalýz algoritmu McEliece
a jeho variant. Zvažte metody zabývající se zkrácením velikosti klíčů.
Implementujte šifrovací a dešifrovací algoritmy a změřte jejich výpočetní
časovou a prostorovou náročnost v závislosti na velikosti klíče.

## Obsah

Stručný obsah repozitáře:

### refs

Ve větvi `master` jsou finální soubory týkající se diplomové práce na téma
*McEliece*. Obsahuje implementaci v softwaru *Mathematica* a text práce ve
formátu *LaTeX*.

Ve větvi `mereni` jsou upravené zdrojové kódy balíků, *Mathematica* "skript"
pro provedení měření a naměřená data.

Tag `odevzdani` označuje orazítkovaný "release" 10.5.2016.

### Adresářová struktura

```
.
+ implementace        implementace McEliece ve Wolfram Mathamatica
|  + grafy              výstupy z měření
|  + src                zdrojové kódy balíků (knihoven)
|  \ *.nb               ukázky a příklady použití
+ materialy           pomocné soubory a skripty
\ text                  text práce
```

