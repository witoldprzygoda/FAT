# Rozkład M(e⁺e⁻) rozpadu η → π⁺π⁻e⁺e⁻ w symulacji SMASH i poprawki w postaci wag

HADES, p+p przy T = 4,53 GeV, kanał π⁺π⁻e⁺e⁻. Notatka z 12 listopada 2026 r. dla
grupy analizy. Gałąź `pp45_pippimepem_betanew_sim`; zmiany nie są jeszcze zapisane
w repozytorium. Liczby z symulacji policzono z regułami wag w wersji 2 (13 listopada
2026 r.), w których stosunki rozgałęzień rozpadów π⁰ i η wykonywanych przez HGeant
pochodzą wprost z konfiguracji HGeant tej produkcji (karta `ETA`, punkt 5 rozdziału 8).

**Oznaczenia.**

- M4 = M(π⁺π⁻e⁺e⁻), OA = kąt otwarcia pary e⁺e⁻.
- **Garb**: η → π⁺π⁻π⁰, para e⁺e⁻ z π⁰ (rozpad Dalitza albo konwersja fotonu),
  jeden foton nie jest rejestrowany. **Pik**: η → π⁺π⁻e⁺e⁻ (rozpad bezpośredni)
  oraz η → π⁺π⁻γ z konwersją fotonu.
- Okno garbu: M4 = 0,455–0,5225 GeV/c²; okno piku: M4 = 0,5336–0,5542 GeV/c².
- Selekcje: `pre` = isBest ≥ 0, PT3, vz > −160 mm, `cor_ok` = 1; E1 = `pre` +
  oa_lab < 50°, kąt w układzie spoczynkowym η > 140°, M(π⁺π⁻) < 0,46 GeV/c²,
  OA > 2° (selekcja z prezentacji); E2 = E1 bez cięcia OA > 2°.
- Błąd liczby zdarzeń piku w danych to błąd statystyczny złożony kwadratowo
  z rozrzutem 48 wariantów modelu dopasowania piku (metoda P_main z badania danych).
- Odchylenie = (symulacja − dane)/σ.
- Każdy stosunek dane/symulacja w skali bezwzględnej zależy od świetlności
  scałkowanej danych L_data = 6080,5 nb⁻¹, która jest **NIEPOTWIERDZONA**: zmienia się
  jak 6080,5/L_data. Takie liczby nie są wynikami fizycznymi.

## 1. Streszczenie

W zestawie η symulacji SMASH rozpad η → π⁺π⁻e⁺e⁻ wygenerowano w PLUTO jako zwykły
rozpad na cztery cząstki, bez dileptonu (fotonu wirtualnego). Dodatek PLUTO z modelami
rozpadów η (`eta_physics`) jest wyłączony celowo, bo narzucałby rozkłady z eksperymentu
DISTO. Dlatego rozkład masy niezmienniczej pary e⁺e⁻ jest niezgodny z QED
już w generatorze: 32 % rozpadów ma M(e⁺e⁻) < 1,52 MeV/c², a tylko 4,7 % ma
M(e⁺e⁻) > 20 MeV/c² (model fizyczny: 6,8 % i 31 %). Akceptancja zmierzona na parach
z rozpadów Dalitza zmienia się z M(e⁺e⁻) łagodnie i tej różnicy nie tłumaczy.
Skutek dla analizy: symulowany pik η ma za dużo par o małym kącie otwarcia.
Wydajność cięcia OA > 2° wynosi 0,590 ± 0,013 wobec 0,778 ± 0,033 w danych, a udział
piku przy M(e⁺e⁻) > 20 MeV/c² 0,046 wobec 0,28 ± 0,07. Bez nowej produkcji
poprawiono to wagami: rozpad bezpośredni przeważono do modelu Krolla–Wady z falą P
pionów, czynnikami postaci i kątami rozpadu, a rozpady η z HGeant (przestrzeń fazowa)
do wykresu Dalitza KLOE-2 (η → π⁺π⁻π⁰) i do czynnika postaci pionu (η → π⁺π⁻γ). Po
poprawkach udział piku przy dużych M(e⁺e⁻) zgadza się z danymi (0,219), rozbieżność
dla cięcia OA > 2° maleje z −5,3σ do −2,3σ, a w garbie poprawia się tylko wydajność
cięcia M(π⁺π⁻) > 0,34 GeV/c². Producentów symulacji prosimy o zdefiniowanie kanału
z modelem pary leptonów i o ponowne wygenerowanie zestawu η (rozdział 8).

## 2. Zestawy symulacji i wagi

Trzy zestawy to konwersje **tych samych** 10¹⁰ zdarzeń SMASH-3.0 (p+p,
E_kin = 4,5 GeV, parametr zderzenia losowany równomiernie w kole o promieniu
1,2 fm). Jedno zdarzenie odpowiada polu π b_max²/N, więc świetlność scałkowana
jednego zestawu wynosi L_sim = 221,05 nb⁻¹. Zestawy różnią się tylko tym, co
konwerter (`convert-smash`) zrobił z niektórymi cząstkami przed HGeant
(opis producentów).

| zestaw | konwerter | co zmieniono | waga wiersza |
|---|---|---|---|
| standardowy (std) | bez PLUTO | nic; π⁰ i η rozpadają się w HGeant | 1, z poprawkami stosunków rozgałęzień HGeant do PDG (np. para z rozpadu Dalitza η: 0,114; piony z η: 1,047) |
| leptonowy (lep) | PLUTO, `database_decays_leptons.dat` | każdy π⁰, η, ρ⁰, ω, φ wymuszony w kanał z parą e⁺e⁻ | (liczba aktywnych kanałów) × BR |
| η (eta) | PLUTO, `database_decays_eta.dat` | każdy η wymuszony: 50 % η → π⁺π⁻e⁺e⁻, 50 % η → X17 π⁺π⁻ (X17 → e⁺e⁻) | produkty mają w = 2; wiersz sygnału W = 2·BR = 5,36·10⁻⁴ przy BR = 2,68·10⁻⁴ (PDG) |

- **Przynależność wierszy.** Każdy rodzaj wiersza bierze się z jednego zestawu
  (podział A, `sim_own_A == 1`). Z zestawu η bierze się tylko wiersze ze śladem
  z wymuszonego rozpadu η. Reszta zestawu η to zestaw standardowy bez naturalnych
  rozpadów η. Uwaga: wiersze zestawu η, które do niego nie należą, mają
  `sim_w` ≈ 1, a nie 0. Zsumowanie ich bez warunku `sim_own_A` dodałoby około
  790 nb, więcej niż całe symulowane widmo (701,8 nb w `pre`).
- **Zgodność normalizacji.** Wiersze bez produktów η w zestawie η i w zestawie
  standardowym dają stosunek 1,0056 ± 0,0036. Wspólna wartość L_sim jest więc
  uzasadniona z dokładnością do 1 %.
- **X17** (klasa 6) służy tylko za kształt przy ilustracyjnym BR = 10⁻⁵ i nigdy nie
  jest dodawany do sumy.
- **Widmo symulowane** to σ = Σ w / L_sim po wierszach std i lep (podział A)
  oraz należących wierszach zestawu η. Sygnałem kanału są wiersze, których
  cząstki według kodów generatora pochodzą dokładnie z tego kanału. Do tego
  dochodzi tło fizyczne. W symulacji nie ma tła kombinatorycznego.
- **Dane**: (N − CB)/L_data, CB = 2√(N₊₊N₋₋) w przedziałach 10 MeV/c².
- **Klasy czwórki** (`src/sim_truth.h`): 2 garb z rozpadem Dalitza π⁰, 3 garb
  z konwersją fotonu z π⁰, 4 pik z konwersją (η → π⁺π⁻γ, zestaw standardowy),
  5 pik η → π⁺π⁻e⁺e⁻ (zestaw η), 6 X17, 7 para e⁺e⁻ z jednego źródła z innymi
  pionami, 8 leptony z różnych źródeł, 0 i 1 błędnie zidentyfikowana cząstka albo
  rozpad pionu.
- **Wagi fizyczne** z rozdziału 5 mnożą `sim_w` danego wiersza. Włącza je opcja
  `SimOpt::phys` (domyślnie wyłączona).

## 3. Spostrzeżenie: rozkład M(e⁺e⁻) rozpadu η → π⁺π⁻e⁺e⁻ w generatorze

W ntuplach są tylko wiersze, które przeszły rekonstrukcję i selekcję. Rozkład
zaakceptowanych wierszy to więc „generator × akceptancja”. Wniosek, że winny jest
generator, opiera się na trzech krokach.

### 3.1. Krok 1: zaakceptowane wiersze a QED

Dla zaakceptowanych wierszy rozpadu bezpośredniego (klasa 5, `pre`) prawdziwą
masę M(e⁺e⁻) liczy się z pędów GEANT obu leptonów. Kształt odniesienia to czynnik
Krolla–Wady razy widmo η → π⁺π⁻γ* z falą P pionów, bez czynników postaci
(`tables/eta_peak_mee_true.txt`; w przedziałach poniżej 25 MeV/c² różni się od
modelu z rozdziału 5 o najwyżej 3 %). Udziały w przedziałach logarytmicznych:

| prawdziwa M(e⁺e⁻) [MeV/c²] | odniesienie QED | zaakceptowane | stosunek |
|---|---:|---:|---:|
| 1,02–1,52 | 0,0665 | 0,3484 | 5,24 |
| 1,52–2,27 | 0,0906 | 0,1828 | 2,02 |
| 2,27–3,37 | 0,0939 | 0,1229 | 1,31 |
| 3,37–5,02 | 0,0945 | 0,0927 | 0,98 |
| 5,02–7,47 | 0,0944 | 0,0723 | 0,77 |
| 7,47–11,13 | 0,0940 | 0,0551 | 0,59 |
| 11,13–16,57 | 0,0931 | 0,0390 | 0,42 |
| 16,57–24,66 | 0,0910 | 0,0276 | 0,30 |
| 24,66–36,71 | 0,0867 | 0,0189 | 0,22 |
| 36,71–54,66 | 0,0783 | 0,0134 | 0,17 |
| 54,66–81,37 | 0,0629 | 0,0097 | 0,15 |
| 81,37–121,14 | 0,0393 | 0,0081 | 0,21 |
| 121,14–180,35 | 0,0138 | 0,0067 | 0,48 |
| 180,35–268,50 | 0,0010 | 0,0026 | 2,59 |
| udział powyżej 24,66 MeV/c² | 0,282 | 0,059 | |

- **Największa gęstość jest tuż nad progiem.** Gęstość zaakceptowanych par na
  MeV/c² jest największa w pierwszym przedziale nad progiem 2m_e (1,022–1,175
  MeV/c²): 1,00 wobec 0,10 w QED. Według QED gęstość dąży do zera na progu.
  Akceptancja nie może wytworzyć takiego maksimum.
- **Spadek jak 1/M².** Między 1,5 a 10 MeV/c² gęstość spada jak M do potęgi −1,73.
  Prawo 1/M² dałoby −2, QED daje około −1 (z dala od progu).
- **Najwyższe przedziały.** Wzrost stosunku przy największych masach pochodzi
  z akceptancji, która dla dużych mas rośnie (rozdział 4.3).

### 3.2. Krok 2: akceptancja zmierzona na parach z rozpadów Dalitza

Generatory rozpadów Dalitza π⁰ i η (HGeant, PLUTO) dają kształt QED. Ten sam
stosunek „zaakceptowane / QED” pokazuje więc dla nich samą akceptancję. Wiersze
z plików pośrednich `sim_compact*.root` (selekcja FAT i vz > −160 mm), bez wag; stosunek unormowano do 1
w przedziale 3,37–5,02 MeV/c² (`tables/rw_dalitz_acceptance.txt`, makro
`sim_rw_dalitz_acc.C`).

| prawdziwa M(e⁺e⁻) [MeV/c²] | Dalitz π⁰, std (HGeant) | Dalitz π⁰, lep (PLUTO) | Dalitz η, lep (PLUTO) | η → π⁺π⁻e⁺e⁻, zestaw η (PLUTO) |
|---|---:|---:|---:|---:|
| 1,02–1,52 | 1,034 | 1,133 | 1,062 | **5,325** |
| 1,52–2,27 | 0,886 | 0,919 | 0,905 | **2,069** |
| 2,27–3,37 | 0,884 | 0,896 | 0,887 | 1,322 |
| 3,37–5,02 | 1 | 1 | 1 | 1 |
| 5,02–7,47 | 1,164 | 1,172 | 1,176 | 0,779 |
| 7,47–11,13 | 1,335 | 1,343 | 1,369 | 0,597 |
| 11,13–16,57 | 1,416 | 1,459 | 1,448 | 0,426 |
| 16,57–24,66 | 1,496 | 1,524 | 1,497 | 0,309 |
| 24,66–36,71 | 1,532 | 1,589 | 1,500 | 0,222 |
| 36,71–54,66 | 1,617 | 1,749 | 1,477 | **0,174** |

- **Pary z rozpadów Dalitza:** trzy próbki z dwóch generatorów dają ten sam łagodny
  przebieg, 0,88–1,75 w zakresie 1–55 MeV/c², w pobliżu progu około 1. Akceptancja
  rośnie z masą, bo pary o większym kącie otwarcia są lepiej rekonstruowane.
- **Rozpad η → π⁺π⁻e⁺e⁻:** przy progu stosunek jest 5 razy większy, przy
  50 MeV/c² 6 razy mniejszy. Zmienia się o czynnik około 30, i to w kierunku
  przeciwnym niż akceptancja.
- **Wniosek:** ten sam detektor i ta sama selekcja nie dają takiej różnicy.
  Pozostaje generator, co potwierdza krok 3.

### 3.3. Krok 3: poziom generatora

Na maszynie produkcyjnej odczytano pliki zdarzeń po konwersji, czyli dane
wejściowe obecnych DST (`/lustre/hades/dstsim/feb22/p4500p/smash/convert/smash_eta/1_old/*.evt`).
Najpierw sprawdzono 40 plików (47 394 rozpady η → π⁺π⁻e⁺e⁻), potem wydobyto próbkę
z 1000 plików (1 181 801 rozpadów; `eta_pipiee_gen_kin_1000files.txt.gz`, rozdział
9). Model docelowy V2 = V3 opisano w rozdziale 5.

| M(e⁺e⁻) [MeV/c²] | generator, 40 plików | generator, 1000 plików | model docelowy (V2 = V3) |
|---|---:|---:|---:|
| 1,02–1,52 | 0,270 ¹ | 0,3203 ¹ | 0,0685 |
| 1,52–2,27 | 0,229 | 0,2289 | 0,0928 |
| 2,27–3,37 | 0,148 | 0,1498 | 0,0964 |
| 3,37–5,02 | 0,098 | 0,1004 | 0,0969 |
| 5,02–7,47 | 0,068 | 0,0670 | 0,0968 |
| 7,47–11,13 | 0,047 | 0,0450 | 0,0962 |
| 11,13–16,57 | 0,031 | 0,0306 | 0,0949 |
| 16,57–24,66 | 0,021 | 0,0205 | 0,0921 |
| 24,66–36,71 | 0,014 | 0,0142 | 0,0865 |
| 36,71–54,66 | 0,011 | 0,0097 | 0,0758 |
| 54,66–81,37 | | 0,0065 | 0,0580 |
| 81,37–121,14 | | 0,0042 | 0,0336 |
| 121,14–180,35 | | 0,0023 | 0,0107 |
| 180,35–268,50 | | 0,0006 | 0,0007 |
| **udział powyżej 20 MeV/c²** | **0,049** | **0,0474** | **0,3133** |

¹ W pliku z próbką 4,94 % wartości M(e⁺e⁻) leży poniżej progu 2m_e, a 3860
wartości to dokładnie 0 (punkt 3.5). W kolumnie z 40 plików te wartości nie
trafiają do żadnego przedziału. W kolumnie z 1000 plików liczy się je do pierwszego
przedziału: 0,3203 − 0,0494 ≈ 0,270, więc obie próbki są zgodne.

- **Kształt generatora.** Rozkład w generatorze maleje mniej więcej jak 1/M².
  W QED rozkład na przedział logarytmiczny jest prawie płaski. Połowa par
  z generatora ma masę poniżej 2,3 MeV/c².
- **Kontrola X17.** Pary z X17 mają w generatorze M(e⁺e⁻) = 17,000 MeV/c² przy
  rozrzucie 0,014 MeV/c², więc ten kanał jest w porządku.

### 3.4. Konfiguracja konwertera

Pliki konwertera odczytano 12 listopada 2026 r. na maszynie produkcyjnej, w katalogu
`/lustre/hades/user/hadesdst/sub/feb22/sim/p4500p/convert_smash/smashtrunk/`.

`database_decays_eta.dat` (wiersze kanałów):

```
type=eta_pi-_pi+           ;channel=eta:e+,e-,pi+,pi-                  ;active=1;
type=eta_X17_pi-_pi+       ;channel=eta:X17[e+,e-],pi+,pi-             ;active=1;
type=eta_pi0dalitz_pi-_pi+ ;channel=eta:pi0[dilepton[e+,e-],g],pi+,pi- ;active=0;
type=eta_pi0_pi-_pi+       ;channel=eta:pi0[g,g],pi+,pi-               ;active=0;
type=eta_g_pi-_pi+         ;channel=eta:g,pi+,pi-                      ;active=0;
```

`pluto_decay.h` (numery wierszy pliku):

```
170:	Ssiz_t i = channel.First(":");       // format examples:  "eta:e+,e-,..."  or  'pi0:dilepton[e+e-],g'
371:	// type=pi0dalitz ;channel=pi0:dilepton[e+,e-],g      ;active=1;
686:	makeDistributionManager()->Disable("eta_physics"); // according to Marten (Iza) , disto physics do messup
693:	makeDistributionManager()->Exec("dalitz_mod: krivoruchenko");
701:	makeStaticData()->AddDecay("X17 --> e+ e-", "X17", "e+, e-", 1);
702:	makeStaticData()->AddDecay("eta --> pi+ pi- X17", "eta", "pi+, pi-, X17", 0.0001);
705:	makeStaticData()->AddDecay("eta --> e+ + e- + pi+ pi-", "eta", "e+,e-,pi+, pi-", 0.0001);
```

- **Aktywny kanał** `eta:e+,e-,pi+,pi-` to zwykły rozpad na cztery cząstki.
  Para e⁺e⁻ nie jest w nim oznaczona jako `dilepton[...]`, choć ta sama baza używa
  tego zapisu np. dla rozpadu Dalitza π⁰ (wiersz 371) i w nieaktywnym kanale
  `eta:pi0[dilepton[e+,e-],g],pi+,pi-`.
- **Definicja rozpadu.** `pluto_decay.h` dodaje rozpad funkcją
  `AddDecay("eta --> e+ + e- + pi+ pi-", ...)` jako ogólny kanał czterociałowy.
  Stosunek rozgałęzień 0,0001 nie ma znaczenia, bo rozpad jest wymuszony.
- **Wyłączony dodatek.** Dodatek PLUTO z modelami rozpadów η (`eta_physics`) jest
  jawnie wyłączony, z komentarzem „according to Marten (Iza), disto physics do
  messup”. Wyłączenie jest celowe: dodatek dawałby rozkłady z eksperymentu DISTO,
  których nie chcemy (informacja od W. Przygody, 13 listopada 2026).

### 3.5. Co dokładnie robi generator (próbka 1 181 801 rozpadów)

- **Rozkład mas.** Gęstość generatora dobrze opisuje postać iloczynowa
  G(q, M) = g(q) · h(M) · λ^{1/2}(m_η², M², q²), gdzie q = M(e⁺e⁻),
  M = M(π⁺π⁻), m_η = 0,548 GeV/c² (wartość PLUTO), a g i h dopasowano metodą
  największej wiarygodności. χ² wynosi 1049 dla 1004 komórek (przybliżona liczba
  stopni swobody 934, p = 0,005). Niezależna, nieparametryczna ocena gęstości daje
  te same wagi w granicach statystyki i zmienia główne liczby o mniej niż 0,3 %.
  PLUTO losowało więc obie masy niezależnie i ważyło je dwuciałową przestrzenią
  fazową. Dlaczego g(q) ma kształt bliski 1/q², nie jest ustalone (NIESPRAWDZONE).
- **Kąty rozpadu** są niezależne od mas i od siebie (korelacja
  0,0002 ± 0,0009).
  - Kąt e⁺ (w układzie spoczynkowym pary e⁺e⁻ względem kierunku pary w układzie η)
    ma rozkład płaski: dopasowanie 1 + b cos²θ daje b = 0,006 ± 0,004 dla
    M(e⁺e⁻) > 2 MeV/c². QED daje 1 + cos²θ (z dala od progu).
  - Kąt π⁺ (w układzie π⁺π⁻) ma już rozkład sin²θ (fala P):
    b = −0,9994 ± 0,0010.
- **Ograniczona dokładność zapisu.** 4,94 % wartości M(e⁺e⁻) w pliku z próbką leży
  poniżej 2m_e, a 3860 to dokładnie 0. Wśród 572 651 symulowanych wierszy klasy 5
  żaden nie ma M(e⁺e⁻) poniżej 2m_e. Rozmycie dotyczy więc zapisanych liczb
  w wydobytej próbce (albo energii w plikach `.evt`), a nie pędów użytych przez
  HGeant. Przyczyna jest NIESPRAWDZONA. W wagach zastosowano jeden przedział
  progowy 2m_e–1,552 MeV/c² dla G i dla modelu. Oszacowany wpływ na sumy rozpadu
  bezpośredniego jest poniżej 1 %.

## 4. Skutki dla analizy (symulacja bez poprawek)

### 4.1. Wydajność cięcia OA > 2°

- **Skąd bierze się zależność od M(e⁺e⁻).** Kąt otwarcia rośnie z M(e⁺e⁻).
  W rozpadzie bezpośrednim cięcie OA > 2° przechodzi tak:

  | prawdziwa M(e⁺e⁻) | udział przechodzący OA > 2° |
  |---|---:|
  | przy progu | 0,525 |
  | 8 MeV/c² | 0,70 |
  | 17 MeV/c² | 0,97 |
  | powyżej 38 MeV/c² | 1,00 |

- **Wydajność OA > 2° w `pre`, wszystkie M4:**

  | | wydajność |
  |---|---:|
  | rozpad bezpośredni | 0,614 |
  | konwersja | 0,564 ± 0,028 |
  | cały symulowany pik | 0,590 ± 0,013 |
  | dane | 0,778 ± 0,033 |

  Odchylenie wynosi −5,3.
- **Cięcie OA < 4°:** 0,852 ± 0,008 w symulacji wobec 0,586 ± 0,056 w danych,
  odchylenie +4,7.

### 4.2. Pik w przedziałach OA i M(e⁺e⁻)

Udziały zawartości piku w selekcji `pre` (dane: P_main w każdym przedziale):

| przedział | dane | symulacja |
|---|---:|---:|
| OA < 2° | 0,220 ± 0,073 | 0,410 |
| OA > 10° | 0,201 ± 0,074 | 0,026 |
| OA > 6° | 0,278 ± 0,079 | 0,053 |
| M(e⁺e⁻) > 20 MeV/c² | 0,280 ± 0,071 | 0,046 |

- **Jeden wspólny czynnik** nie opisuje przedziałów: χ²/ndf = 9,5/4 (p = 0,05)
  dla OA i 10,8/4 (p = 0,03) dla M(e⁺e⁻).
- **Osobne czynniki** dla konwersji i rozpadu bezpośredniego wymagają ujemnego
  czynnika konwersji: −1,07 ± 0,50 (OA) i −1,53 ± 0,56 (M(e⁺e⁻)). Żadna fizyczna
  mieszanina dwóch symulowanych kształtów nie opisuje danych.
- **Moc testu.** Błędy danych to 30–40 % na przedział, więc testy na przedziałach
  mają niewielką moc. Najwyraźniejszym sprawdzianem jest udział piku przy dużych
  M(e⁺e⁻) i dużych OA.

### 4.3. Normalizacja piku

- **Akceptancja zależy od M(e⁺e⁻).** Względna akceptancja rozpadu bezpośredniego
  rośnie z M(e⁺e⁻): od około 0,8 przy 2 MeV/c² do 1,35 przy 50 MeV/c² i 2,3 przy
  120 MeV/c². Liczba zaakceptowanych rozpadów zależy więc od kształtu
  w generatorze. Na jeden wygenerowany rozpad X17 (M(e⁺e⁻) = 17 MeV/c²) jest
  akceptowany 1,379 razy częściej niż rozpad η → π⁺π⁻e⁺e⁻ z generatora.
- **Udział rozpadu bezpośredniego.** Bez poprawek rozpad bezpośredni daje
  w oknie piku 1,13 (`pre`) i 1,30 (E1) raza tyle co konwersja. To nie jest
  przewidywanie fizyczne.
- **Zależność od metody.** Stosunek dane/symulacja dla piku zależy od metody
  wyznaczenia: 0,29–1,36 dla metod A–E (`tables/eta_summary_peak_methods.txt`).
  Normalizacja piku jest znana tylko z dokładnością do czynnika około 2. Pewne jest
  tylko to, że pik w danych nie jest większy niż w symulacji.
- **Wcześniejsza ilustracja.** Przeważenie zaakceptowanych wierszy do kształtu QED
  przy założeniu płaskiej akceptancji dawało wydajność OA > 2° równą 0,775. Była to
  tylko ilustracja kształtów; jej wyników nie należy cytować. Poprawne przeważenie
  na poziomie generatora opisuje rozdział 5.

Osobna sprawa tej samej natury: w zestawie standardowym rozpady η → π⁺π⁻π⁰
i η → π⁺π⁻γ wykonuje HGeant według przestrzeni fazowej (punkt 5.7). Wydajność
cięcia M(π⁺π⁻) > 0,34 GeV/c² dla garbu znakowanego fotonem wynosi w symulacji
0,444 ± 0,025, a w danych 0,705 ± 0,017.

## 5. Poprawki fizyczne zastosowane jako wagi

Decyzja: poprawki w postaci wag, bez nowej produkcji. Każdy wiersz dostaje czynnik,
który mnoży jego wagę `sim_w`. Wzory razem ze źródłami są w
`plots/eta_sim/sim_eta_physics.h`. Wagi zależą tylko od prawdziwych (GEANT) mas
i kątów, nie od wielkości zrekonstruowanych.

### 5.1. Rozpad bezpośredni η → π⁺π⁻e⁺e⁻ (klasa 5)

```
w5 = T / G                       q = M(e⁺e⁻), M = M(π⁺π⁻)   (wartości prawdziwe)

T(q, M, θe, θπ) = KW(q) · R(q, M) · A_e(θe; q) · sin²θπ

KW(q)    = (2α / (3π q)) · (1 + 2m_e²/q²) · √(1 − 4m_e²/q²)
R(q, M)  = p*³ · λ^{3/2}(m_η², M², q²) · |F_V(M²)|² · |F_VMD(q²)|²
A_e(θe)  = 1 + cos²θe + (4m_e²/q²) · sin²θe
λ(a,b,c) = a² + b² + c² − 2ab − 2ac − 2bc
p*       = pęd pionu w układzie spoczynkowym π⁺π⁻ = M β_π / 2
F_V      = czynnik postaci pionu, postać Guerrero–Picha (wariant: Gounaris–Sakurai);
           m_ρ = 0,7753 GeV/c², Γ_ρ = 0,1491 GeV
F_VMD    = m_ρ² / (m_ρ² − q²)   (biegun ρ; szerokość 0 poniżej 4m_π², czyli w całym
           zakresie q < m_η − 2m_π = 0,269 GeV/c²)

G(q, M, θe, θπ) = g(q) · h(M) · λ^{1/2}(m_η², M², q²) · (kąt e⁺ płaski) · sin²θπ
```

- **Źródła:**
  - czynnik Krolla–Wady: N. M. Kroll, W. Wada, Phys. Rev. 98 (1955) 1355;
  - F_V: F. Guerrero, A. Pich, Phys. Lett. B 412 (1997) 382, hep-ph/9707347;
  - wariant F_V: G. J. Gounaris, J. J. Sakurai, Phys. Rev. Lett. 21 (1968) 244.
- **Skąd p*³ λ^{3/2}.** Amplituda anomalna η → π⁺π⁻γ* ma postać
  ε_{μναβ} ε^μ p₊^ν p₋^α k^β (fala P pionów, foton poprzeczny), więc
  |A|² ∝ p*² λ sin²θπ. Przestrzeń fazowa daje λ^{1/2} p* dM, stąd
  dΓ/dM ∝ p*³ λ^{3/2}. Dla q → 0 to przechodzi w znane widmo η → π⁺π⁻γ,
  proporcjonalne do p*³ k³ |F_V|² (k: pęd fotonu w układzie η).
- **Kąty.** θe to kąt e⁺ w układzie spoczynkowym pary e⁺e⁻ względem kierunku pary
  w układzie η, a θπ to kąt π⁺ w układzie π⁺π⁻ względem kierunku pary π⁺π⁻.
  - Pełny rozkład kątowy tej samej amplitudy to sin²θπ · [1 − β_e² sin²θe sin²φ],
    gdzie φ to kąt między płaszczyznami rozpadu obu par, a β_e² = 1 − 4m_e²/q².
    A_e jest jego średnią po φ.
  - Generator ma już rozkład sin²θπ, więc czynnik pionowy wagi wynosi 1.
- **Kąty wierszy.** Prawdziwe kąty każdego symulowanego wiersza policzono z pędów
  GEANT w `scan_nt`. Znaleziono je dla wszystkich 572 651 wierszy klasy 5;
  0 przypadków nie dało się dopasować.
- **Porównanie z analizą X17.** W prezentacji analizy X17 tego kanału (K. Proscinski)
  zastosowano poprawki tego samego rodzaju: czynnik VMD dla η → π⁺π⁻e⁺e⁻ i czynnik
  postaci pionu dla η → π⁺π⁻γ.

### 5.2. Garb: wykres Dalitza η → π⁺π⁻π⁰ (klasy 2 i 3)

```
f(X, Y) = 1 + aY + bY² + dX² + fY³ + gX²Y
a = −1,095   b = 0,145   d = 0,081   f = 0,141   g = −0,044         (KLOE-2)
X = √3 (T₊ − T₋) / Q,   Y = 3T₀/Q − 1,   Q = T₊ + T₋ + T₀
T₀ = (m_η² + m_π0² − M²) / (2m_η) − m_π0      (energia kinetyczna π⁰ w układzie η)
X² → ⟨X²⟩_Y = X_max(Y)² / 3,   X_max = 2√3 · p* · q₀* / (m_η Q)
wB = f̄(Y) / ⟨f⟩_PS,   ⟨f⟩_PS = 1,0384
```

- **Źródło:** A. Anastasi i in. (KLOE-2), JHEP 05 (2016) 019, arXiv:1601.06985.
  p* i q₀* to pędy π± i π⁰ w układzie spoczynkowym π⁺π⁻.
- **Dlaczego uśrednianie po X.** Foton z π⁰ nie jest rejestrowany, więc znamy tylko
  M(π⁺π⁻), czyli Y. Przy ustalonym Y przestrzeń fazowa jest płaska w X między
  ±X_max, dlatego X² zastępuje się średnią X_max²/3.
  - Sprawdzenie symulacją przestrzeni fazowej (4·10⁶ rozpadów): ⟨X²⟩ zgadza się
    z X_max²/3 w granicach 0,3 %.
  - ⟨f⟩_PS z pełnym f(X, Y) wynosi 1,03838, z zastąpieniem X² 1,03837.
- **Wariant z analizy X17:** f = 1 + aY, a = −1,095, ⟨f⟩_PS = 0,9900.
- **Zakres wag.** Wagi wierszy klas 2 i 3 w `pre` mają średnią 0,889 i 0,892
  i zakres 0,23–2,01. Powyżej M(π⁺π⁻) = m_η − m_π0 = 0,413 GeV/c² garbu nie ma.

### 5.3. η → π⁺π⁻γ (pik z konwersją, klasa 4)

```
wC = |F_V(M²)|² · p*² · k² / ⟨|F_V|² p*² k²⟩_PS
```

- **Skąd ta postać.** HGeant daje przestrzeń fazową, czyli dΦ/dM ∝ p* k (k: pęd
  fotonu w układzie η). Rozkład fizyczny to p*³ k³ |F_V|², a ich stosunek jest
  powyżej.
- **Stała normalizacyjna:** 5,229·10⁻⁴ (F_V Guerrero–Picha), policzona numerycznie
  i sprawdzona symulacją przestrzeni fazowej.
- **Zakres wag.** Wagi 307 wierszy klasy 4 w `pre` mają średnią 1,014 i zakres
  0,007–1,39.

### 5.4. Pozostałe wiersze

- **Zestaw η, para e⁺e⁻ z wymuszonego rozpadu poza klasą 5** (klasa 7 i inne
  wiersze z klasą pary leptonów L10): tylko waga w M(e⁺e⁻), T_ee/G_ee (rozkłady
  brzegowe T i G). Czynnika kątowego nie ma, także w wariancie V3, bo kątów tych
  wierszy nie da się policzyć.
- **Zestaw η, piony z wymuszonego rozpadu poza klasą 5** (klasa pary pionów P9):
  waga T_ππ/G_ππ.
- **Zestawy std i lep, piony z jednego η poza klasami 2–4** (klasa pary pionów P2;
  1547 wierszy w `pre`): mieszanina p3·wB + (1 − p3)·wC.
  - p3 to udział π⁺π⁻π⁰ przy danym M(π⁺π⁻), liczony ze stosunków rozgałęzień HGeant
    0,2114 (π⁺π⁻π⁰) i 0,0478 (π⁺π⁻γ) oraz z gęstości przestrzeni fazowej.
  - Powyżej m_η − m_π0 możliwy jest tylko π⁺π⁻γ.
- **X17 (klasa 6):** bez wagi.

### 5.5. Normalizacja

- **Liczba rozpadów bezpośrednich bez zmian.** T i G są unormowane do 1 na całej
  przestrzeni fazowej. Średnia w5 po 1 181 801 rozpadach z generatora wynosi przed
  ostatnim przemnożeniem 1,00003 (V3) i 0,99978 (V2), a po nim dokładnie 1.
  Liczba rozpadów i BR = 2,68·10⁻⁴ się nie zmieniają.
- **wB i wC też nie zmieniają liczby rozpadów.** Ich średnia po przestrzeni
  fazowej wynosi 1. Zmienia się tylko kształt, a liczba rozpadów η → π⁺π⁻π⁰
  i η → π⁺π⁻γ zostaje. Zmiany przekrojów czynnych zaakceptowanych wierszy
  (rozdział 6) wynikają z akceptancji.
- **Sprawdzenie na próbce z generatora.** Przeważona próbka zgadza się z T:

  | wariant | M(e⁺e⁻) | M(π⁺π⁻) | kąt e⁺ (M(e⁺e⁻) > 2 MeV/c²) |
  |---|---:|---:|---:|
  | V3 | χ² = 14,3/30 | 28,8/27 | 23,0/20 |
  | V2 | 9,3/30 | 27,1/27 | |

  Wagi V3 na tej próbce mają rms 1,57 i zakres 0,004–16,8; efektywna liczba
  zdarzeń wynosi 0,29 liczby rozpadów.
- **Stosunek do PDG.** Model daje Γ(π⁺π⁻e⁺e⁻)/Γ(π⁺π⁻γ) = 6,33·10⁻³ (V2 = V3),
  wobec (6,35 ± 0,29)·10⁻³ według PDG. To tylko sprawdzenie zgodności normalizacji:
  stosunek zależy głównie od małych M(e⁺e⁻), a F_VMD zmienia go o 0,4 %. Kształtu
  M(e⁺e⁻) ten test nie sprawdza.
- **Niezależna kontrola.** Recenzent zaprogramował wszystkie wzory od nowa
  i otrzymał te same wartości:
  - wB zgadza się do 1,3·10⁻⁴ we wszystkich 2836 wierszach garbu;
  - wC zgadza się do 7·10⁻⁸ we wszystkich 307 wierszach konwersji;
  - wydajność OA > 2° poprawionego piku: 0,6996 ± 0,0143;
  - garb w oknie garbu w E1: 5,73 ± 0,18 nb.

  Recenzent liczył z regułami wag w wersji 1. Z regułami w wersji 2 makro daje
  0,697 ± 0,015 i 5,64 ± 0,17 nb.

### 5.6. Warianty

| wariant | rozpad bezpośredni | uwagi |
|---|---|---|
| V0 | KW · p*³ λ^{3/2} | bez czynników postaci i kątów |
| V1 | V0 · \|F_V\|² | |
| V2 | V1 · \|F_VMD\|² | główny przed recenzją; zostawia płaski kąt e⁺ z generatora |
| **V3 (główny)** | V2 · A_e · sin²θπ | pełny model tej samej amplitudy po uśrednieniu po φ |
| V4 | V3 · człon płaszczyzn R_φ | tylko sprawdzenie |
| GS | V3 z F_V Gounarisa–Sakuraia | |
| 1 + aY | V3, garb z f = 1 + aY | wariant z analizy X17 |

Garb (KLOE-2) i konwersja (F_V Guerrero–Picha) są takie same we wszystkich
wariantach poza „1 + aY” (inny garb) i GS (konwersja też z F_V
Gounarisa–Sakuraia). Wpływ wariantów na główne liczby (`tables/rw_compare_variants.txt`):

- **GS:** najwyżej około 1 % od V3; **V1:** najwyżej około 1 % od V2 (oba bez kątów).
- **V0:** do 14 % od V3 (symulowany pik w `pre`).
- **V2 względem V3:**
  - zaakceptowany rozpad bezpośredni jest w V3 mniejszy o 14 % (V3/V2 = 0,862
    w `pre`). Akceptancja jest największa dla kątów e⁺ w pobliżu 90°, a czynnik
    A_e przesuwa pozytony ku cos θe = ±1;
  - stosunki dane/symulacja dla piku rosną w V3 o około 9 %.
- **V4:** V4/V3 = 0,998–1,000 we wszystkich selekcjach, ⟨cos 2φ⟩ = −0,0020 ± 0,0017
  na zaakceptowanych wierszach. Członu płaszczyzn można nie uwzględniać.
- **1 + aY:** wydajność cięcia M(π⁺π⁻) > 0,34 dla garbu znakowanego fotonem rośnie
  o 5 %; pozostałe liczby zmieniają się o najwyżej 2 %, a χ² dopasowania M1 w E1
  rośnie z 210 do 218.

### 5.7. Co zakłada się o HGeant

Wagi garbu i konwersji zakładają, że HGeant wykonuje rozpady η → π⁺π⁻π⁰
i η → π⁺π⁻γ według przestrzeni fazowej. W lokalnym źródle HGeant 2-7.4
(`/home/przygoda/hydra2/hydra2-git/hgeant2-7.4/geant/gdecay.F`) tylko η → γe⁺e⁻
idzie do GDALET. Każdy inny trzyciałowy rozpad η idzie do GDECA3, opisanego
komentarzem „phase space”:

```
         ELSE IF (IPART.EQ.17) THEN     ! eta
            IF(MXX.EQ.30201) THEN       ! -> gamma e+e-
               CALL GDALET(DMASS,XM(1),XM(2),XM(3),PCM,9996)
            ELSE    ! phase space
               CALL GDECA3(DMASS,XM(1),XM(2),XM(3),PCM)
```

Założenie jest więc sprawdzone w kodzie źródłowym. NIESPRAWDZONE jest to, że
produkcja użyła tej wersji HGeant.

### 5.8. Użycie w kodzie

```cpp
#include "sim_eta_physics.h"            // po SimEngine.h
simphys::Init(); simphys::Install();    // wczytuje G i stałe, podłącza wagę do silnika
SimOpt o; o.phys = simphys::kPhysMain;  // 108: V3 + garb KLOE-2 + konwersja F_V + zasady dla tła
```

- **Wartość domyślna** `o.phys = 0` wyłącza wagi; wtedy wszystkie dotychczasowe
  wyniki są bez zmian.
- **Dawny wariant główny:** `simphys::kPhysV2` = 107.
- **Kąty V3** wczytuje `simphys::LoadAngles()`. Gdy się o tym zapomni, V3 wczyta
  je sam przy pierwszym użyciu.

## 6. Wyniki po poprawkach

Porównanie przed i po wagach liczy `plots/eta_sim/sim_rw_compare.C`. Makro
sprawdza 103 liczby z wcześniejszych tabel; z regułami wag w wersji 2 (po
dostosowaniu wartości wzorcowych 13 listopada 2026) wszystkie mieszczą się w tolerancji.
„Po” oznacza wariant główny V3. Kolumna V2 jest podana tam, gdzie różnica ma
znaczenie.

### 6.1. Wielkości niezależne od skali bezwzględnej (to są wyniki)

W nawiasach podano odchylenie (symulacja − dane)/σ.

| wielkość | selekcja | dane | symulacja przed | po (V3) | po (V2) |
|---|---|---:|---:|---:|---:|
| udział piku z M(e⁺e⁻) > 20 MeV/c² | `pre` | 0,280 ± 0,071 | 0,046 (−3,3) | **0,219 (−0,8)** | 0,235 (−0,6) |
| | `pre` + OA > 2° | 0,358 ± 0,085 | 0,077 (−3,3) | 0,315 (−0,5) | 0,332 (−0,3) |
| | E1 | 0,359 ± 0,095 | 0,085 (−2,9) | 0,330 (−0,3) | 0,347 (−0,1) |
| | E2 | 0,284 ± 0,080 | 0,048 (−3,0) | 0,225 (−0,7) | 0,241 (−0,5) |
| udział piku z OA > 6° | `pre` | 0,278 ± 0,079 | 0,053 (−2,8) | 0,222 (−0,7) | 0,239 (−0,5) |
| | `pre` + OA > 2° | 0,356 ± 0,096 | 0,090 (−2,8) | 0,318 (−0,4) | 0,337 (−0,2) |
| udział piku z OA < 2° | `pre` | 0,220 ± 0,073 | 0,410 | 0,303 | |
| wydajność OA > 2° (pik) | `pre` | 0,778 ± 0,033 | 0,590 ± 0,013 (−5,3) | **0,697 ± 0,015 (−2,3)** | 0,708 ± 0,013 (−2,0) |
| wydajność OA < 4° (pik) | `pre` | 0,586 ± 0,056 | 0,852 ± 0,008 (+4,7) | 0,671 ± 0,011 (+1,5) | 0,651 ± 0,010 (+1,1) |
| wydajność M(π⁺π⁻) > 0,34 (pik) | `pre` | 0,845 ± 0,050 | 0,756 ± 0,013 (−1,7) | 0,816 ± 0,009 (−0,6) | 0,816 ± 0,008 (−0,6) |
| wydajność kąta w układzie η > 160° (pik) | `pre` | 0,918 ± 0,049 | 0,993 ± 0,002 (+1,5) | 0,992 ± 0,002 (+1,5) | |
| wydajność pT4 > 400 MeV/c (pik) | `pre` | 0,690 ± 0,064 | 0,797 ± 0,010 (+1,7) | 0,804 ± 0,011 (+1,8) | 0,801 ± 0,010 |
| wydajność oa_lab < 50° (pik) | `pre` | 0,920 ± 0,035 | 0,901 ± 0,009 (−0,5) | 0,918 ± 0,008 (−0,1) | |
| wydajność „oba piony w dopasowaniu wierzchołka” (pik) | `pre` | 0,940 ± 0,041 | 0,927 ± 0,007 (−0,3) | 0,930 ± 0,007 (−0,2) | |
| jeden czynnik na przedziałach: χ²/ndf, OA / M(e⁺e⁻) | `pre` | | 9,5/4 / 10,8/4 | 1,7/4 / 1,0/4 | 1,3/4 / 0,8/4 |
| | E2 | | 14,8/4 / 11,6/4 | 4,0/4 / 1,1/4 | 3,4/4 / 0,8/4 |
| dwa swobodne czynniki (przedziały OA): konwersja / rozpad bezpośredni | `pre` | | −1,07 ± 0,50 / 1,47 ± 0,47 | 0,17 ± 0,12 / 0,39 ± 0,11 | 0,18 ± 0,12 / 0,33 ± 0,10 |
| garb znakowany fotonem: wydajność M(π⁺π⁻) > 0,34 | `pre` | 0,705 ± 0,017 | 0,444 ± 0,025 (−8,6) | 0,667 ± 0,024 (−1,3); wariant 1 + aY: 0,700 ± 0,023 (−0,2) | |
| kształt M(π⁺π⁻) w oknie garbu przy ustalonym składzie: χ²/ndf | E1 / E2 / `pre` | | 69,4/18 / 109,8/18 / 79,5/19 | 83,9/18 / 124,9/18 / 217,2/19 (gorzej) | |
| położenie piku [GeV/c²] / szerokość [MeV/c²] | `pre` | 0,54296 ± 0,00046 / 5,15 ± 0,55 | 0,54602 ± 0,00019 / 6,32 ± 0,18 | 0,54629 ± 0,00019 / 6,00 ± 0,16 | |

- **Cięcia M(π⁺π⁻).** W danych wydajności cięć M(π⁺π⁻) dla piku są zniekształcone
  przez samą selekcję (wniosek z badania danych).
- **Cięcie kąta w układzie η > 160°.** Dopasowania próbki odrzuconej przez to cięcie
  mają we wszystkich wariantach p < 0,05.
- **Znakowanie fotonem.** W symulacji jest naśladowane warunkiem na foton
  0,119 ≤ M(e⁺e⁻γ) < 0,173 GeV/c² (NIESPRAWDZONE). Wydajność w danych pochodzi
  z η znakowanych w ECAL, z badania danych bez cięcia vz.
- **Liczebność próbek.** Swobodne czynniki pochodzą z małych wzorców:
  - konwersja: 307 wierszy (efektywnie 270);
  - garb: 2836 wierszy (efektywnie 2279).

  Czynnik konwersji jest tylko około 1,5σ od zera. Oba czynniki są dodatnie
  i zgodne, ale nie mierzą podziału piku w danych na konwersję i rozpad
  bezpośredni.

### 6.2. Przekroje czynne w symulacji (Σ w / L_sim, bez poprawki na wydajność)

| składnik | selekcja, okno | przed [nb] | po (V3) [nb] | po/przed |
|---|---|---:|---:|---:|
| rozpad bezpośredni | `pre`, wszystkie M4 | 1,3852 ± 0,0018 | 1,4168 ± 0,0034 | 1,023 (V2: 1,187) |
| rozpad bezpośredni | `pre` + OA > 2°, wszystkie M4 | 0,8501 ± 0,0014 | 1,1495 ± 0,0033 | 1,352 |
| konwersja w piku | `pre`, wszystkie M4 | 1,225 ± 0,070 | 1,242 ± 0,076 | 1,014 |
| garb (klasy 2 + 3) | `pre`, wszystkie M4 | 13,75 ± 0,26 | 12,24 ± 0,26 | 0,890 |
| garb (klasy 2 + 3) | E1, okno garbu | 6,17 ± 0,17 | 5,64 ± 0,17 | 0,913 |
| pik (klasy 4 + 5) | E1, okno piku | 0,980 ± 0,041 | 1,272 ± 0,048 | 1,297 |
| pik (klasy 4 + 5) | E2, okno piku | 1,766 ± 0,057 | 1,881 ± 0,065 | 1,065 |
| suma T1 + T2 | `pre`, wszystkie M4 | 701,8 ± 1,5 | 701,1 ± 1,5 | 0,999 |

Wagi zachowują liczbę rozpadów w generatorze. Wzrost zaakceptowanego rozpadu
bezpośredniego, szczególnie z cięciem OA > 2°, pochodzi więc z akceptancji, która
rośnie z M(e⁺e⁻).

### 6.3. Wielkości bezwzględne (to NIE są wyniki fizyczne)

Wszystkie poniższe liczby zmieniają się jak 6080,5/L_data (L_data NIEPOTWIERDZONA).
Stosunek dla piku zależy dodatkowo od metody: metody A–E różnią się po poprawkach
o czynnik 3,7.

| wielkość | przed | po (V3) |
|---|---:|---:|
| pik, dane/symulacja, metody A–E (`pre`, `pre` + OA > 2°, E1, E2) | 0,29–1,36 | 0,28–1,02 |
| dane/symulacja w oknie piku, E1 / E2 | 1,264 ± 0,029 / 1,155 ± 0,022 | 1,184 ± 0,026 / 1,127 ± 0,022 |
| dane/symulacja w oknie garbu, E1 / E2 | 0,953 ± 0,016 / 0,901 ± 0,013 | 0,984 ± 0,017 / 0,941 ± 0,014 |
| dane/symulacja poniżej garbu (M4 = 0,35–0,455), E1 / E2 | 1,037 ± 0,061 / 0,960 ± 0,048 | 1,294 ± 0,076 / 1,185 ± 0,060 |
| skok na krawędzi garbu (pas 0,56–0,595 / pas 0,4875–0,5225), E1 | 1,466 | 1,437 ± 0,041 |
| dopasowanie wzorców M1 w E1 (0,40–0,75 GeV/c²): α_garb / α_pik / χ² (67 st. swobody) | 0,44 ± 0,07 / 1,04 ± 0,27 / 193 | 0,49 ± 0,08 / 0,73 ± 0,22 / 210 |

### 6.4. Co poprawki naprawiają, a czego nie

**Naprawione:**

- **Udział piku przy dużych M(e⁺e⁻) i dużych OA zgadza się z danymi** (`pre`:
  0,219 wobec 0,280 ± 0,071; przed poprawkami 0,046). To główny wynik.
- **Przedziały OA i M(e⁺e⁻) opisuje jeden wspólny czynnik**, a swobodne czynniki
  konwersji i rozpadu bezpośredniego są dodatnie. Testy na przedziałach mają
  jednak niewielką moc.
- **Rozbieżność dla OA < 4°** maleje z +4,7 do +1,5.
- **Garb:** wydajność cięcia M(π⁺π⁻) > 0,34 dla garbu znakowanego fotonem
  (−8,6 → −1,3).

**Nie naprawione:**

- **OA > 2°:** rozbieżność maleje z −5,3 do −2,3 (V2: −2,0), ale nie znika.
- **Kąt w układzie η > 160° i pT4 > 400 MeV/c:** odchylenia +1,5 i +1,8 bez zmian.
- **Garb:** waga KLOE-2 nie naprawia kształtu.
  - Przy ustalonym składzie symulacji kształt M(π⁺π⁻) w oknie garbu opisuje się
    gorzej.
  - Udział powyżej 0,34 GeV/c² jest teraz za duży: w E1 dane 0,556 ± 0,002,
    symulacja 0,510 → 0,591.
  - Poniżej garbu otwiera się niedobór symulacji.
  - Poprawiony garb razem z tłem opisuje M(π⁺π⁻) i M4 tylko wtedy, gdy garb
    przemnoży się przez czynnik zależny od badanej wielkości: 0,28–0,40
    z dopasowania M(π⁺π⁻), 0,44–0,49 z dopasowania M4 (słabe dopasowania; to nie
    są liczby zdarzeń).
  - Skok na krawędzi garbu zostaje (E1: 1,47 → 1,44).
  - Nie wykazano, że problem garbu dotyczy tylko normalizacji.
- **Dopasowania wzorców w M4 nie są lepsze** (E1: χ² = 193 → 210 na 67 stopni
  swobody), a normalizacja piku nadal zależy od metody.

### 6.5. Rysunki

Rysunki nie są w repozytorium; odtwarzają je makra z rozdziału 9.

- [eta_peak_mee_true.png](../plots/eta_sim/figs/eta_peak_mee_true.png): prawdziwa
  M(e⁺e⁻) zaakceptowanego rozpadu bezpośredniego wobec kształtu QED (krok 1).
- [rw_weights_generator.png](../plots/eta_sim/figs/rw_weights_generator.png): M(e⁺e⁻)
  i M(π⁺π⁻) w generatorze wobec modeli V0–V2, stosunek model/generator, test postaci
  iloczynowej i rozkłady kątów e⁺ i π⁺ (krok 3).
- [rw_weights_acceptance.png](../plots/eta_sim/figs/rw_weights_acceptance.png):
  względna akceptancja rozpadu bezpośredniego w funkcji prawdziwej M(e⁺e⁻) i kształty
  M(e⁺e⁻) przed i po wadze.
- [rw_weights_model.png](../plots/eta_sim/figs/rw_weights_model.png): czynniki postaci
  F_V i F_VMD, udziały M(e⁺e⁻) i sprawdzenie stosunku do PDG.
- [rw_weights_closure.png](../plots/eta_sim/figs/rw_weights_closure.png): próbka
  z generatora po przeważeniu wobec modelu T w M(e⁺e⁻), M(π⁺π⁻) i obu kątach.
- [rw_weights_bump.png](../plots/eta_sim/figs/rw_weights_bump.png): waga garbu wB
  w funkcji M(π⁺π⁻) (KLOE-2, 1 + aY), sprawdzenie ⟨X²⟩, waga konwersji wC
  i udział p3.
- [rw_weights_classes.png](../plots/eta_sim/figs/rw_weights_classes.png): rozkłady
  wag w klasach i prawdziwe M(π⁺π⁻) garbu i konwersji przed i po.
- [rw_weights_xsec.png](../plots/eta_sim/figs/rw_weights_xsec.png): stosunki przekrojów
  czynnych po/przed dla składników, selekcji, okien i wariantów.
- [rw_weights_peak.png](../plots/eta_sim/figs/rw_weights_peak.png): zaakceptowany rozpad
  bezpośredni w M(e⁺e⁻) i OA dla wariantów V0–V3 oraz przedziały piku wobec danych.
- [rw_compare_peak_classes.png](../plots/eta_sim/figs/rw_compare_peak_classes.png):
  udziały piku w przedziałach OA i M(e⁺e⁻) w `pre` i E2, dane wobec symulacji przed
  i po, oraz czynniki dopasowań.
- [rw_compare_peak_eff.png](../plots/eta_sim/figs/rw_compare_peak_eff.png): siedem
  wydajności pojedynczych cięć dla piku, dane wobec symulacji przed i po,
  z odchyleniami.
- [rw_compare_peak_position.png](../plots/eta_sim/figs/rw_compare_peak_position.png):
  kształt, położenie i szerokość symulowanego piku wobec danych.
- [rw_compare_bump.png](../plots/eta_sim/figs/rw_compare_bump.png): garb przed i po
  wadze KLOE-2 (i 1 + aY), M(π⁺π⁻) w oknie garbu wobec danych oraz wydajności cięć
  dla garbu znakowanego fotonem.
- [rw_compare_m4_E12.png](../plots/eta_sim/figs/rw_compare_m4_E12.png): dane N, CB
  i N − CB wobec poprawionej symulacji w E1 i E2, ze stosunkami przed i po.
- [rw_compare_m4_E34.png](../plots/eta_sim/figs/rw_compare_m4_E34.png): to samo dla
  selekcji E3 i E4, bez cięcia OA > 2° i z nim.
- [rw_compare_regions.png](../plots/eta_sim/figs/rw_compare_regions.png): dane/symulacja
  w obszarach M4 i skok na krawędzi garbu przed i po.
- [rw_compare_fits.png](../plots/eta_sim/figs/rw_compare_fits.png): czynniki dopasowań
  wzorców M1 i M3 i χ²/ndf przed i po.
- [rw_compare_methods.png](../plots/eta_sim/figs/rw_compare_methods.png): stosunek
  dane/symulacja dla piku według metod A–E przed i po.
- [rw_compare_variants.png](../plots/eta_sim/figs/rw_compare_variants.png): główne
  liczby dla wszystkich wariantów wag.

Trzy najważniejsze:

![Pik η w przedziałach OA i M(e⁺e⁻) przed i po poprawkach](../plots/eta_sim/figs/rw_compare_peak_classes.png)

![Wydajności cięć dla piku η przed i po poprawkach](../plots/eta_sim/figs/rw_compare_peak_eff.png)

![Rozkłady w generatorze wobec modelu docelowego](../plots/eta_sim/figs/rw_weights_generator.png)

## 7. Ograniczenia i rzeczy niesprawdzone

- **L_data = 6080,5 nb⁻¹ jest NIEPOTWIERDZONA.** Każdy bezwzględny stosunek
  dane/symulacja zmienia się jak 6080,5/L_data. Skala bezwzględna jest znana
  z dokładnością około 15–20 %.
- **HGeant liczy trzyciałowe rozpady η według przestrzeni fazowej.** Sprawdzono to
  w lokalnym źródle HGeant 2-7.4. Nie sprawdzono (NIESPRAWDZONE), czy produkcja
  użyła tej wersji.
- **Dokładność zapisu w próbce z generatora.** 4,9 % wartości leży poniżej 2m_e.
  Przyczyna jest NIESPRAWDZONA, a wpływ oszacowano na mniej niż 1 %. Nie zrobiono
  wariantu z przesuniętą granicą przedziału progowego (1,4–2,0 MeV/c²).
- **Kąt φ między płaszczyznami.** Założono, że w generatorze jest płaski; plik
  z próbką nie zawiera φ (NIESPRAWDZONE). Człon płaszczyzn zmienia wyniki
  o mniej niż 0,2 %.
- **Klasa 7 i inne wiersze zestawu η z parą z wymuszonego rozpadu** dostają tylko
  wagę w M(e⁺e⁻), bez czynnika kątowego, także w V3. Oszacowanie wpływu: najwyżej
  około 0,2 % sumy T1 + T2 w oknach (NIESPRAWDZONE).
- **Mieszanina dla wierszy z klasą pionów P2.**
  - Stosuje stosunki rozgałęzień HGeant, a zapisana waga zawiera uśrednioną
    poprawkę 1,047 zamiast osobnych 1,084 (π⁺π⁻π⁰) i 0,882 (π⁺π⁻γ).
  - Wpływ to około 0,15 % w oknie garbu w E1. Oszacowano go z regułami wag
    w wersji 1 (poprawki 1,056, 1,105 i 0,849); z regułami w wersji 2 go nie
    przeliczono (NIESPRAWDZONE).
  - Sama reguła mieszaniny jest NIESPRAWDZONA.
- **Brakujący czynnik w modelu M(π⁺π⁻).** Model M(π⁺π⁻) dla η → π⁺π⁻γ(*) to |F_V|²
  razy przestrzeń fazowa fali P. Pomiary wymagają dodatkowego czynnika (1 + α s)²
  z α około 1,3–1,9 GeV⁻² (wartości recenzent podał z pamięci, NIESPRAWDZONE).
  Nie wprowadzono go. Według recenzenta zmieniłby:
  - udział konwersji z M(π⁺π⁻) > 0,40 GeV/c² z 0,542 na 0,577;
  - Γ(π⁺π⁻e⁺e⁻)/Γ(π⁺π⁻γ) z 6,331·10⁻³ na 6,275·10⁻³.
- **Postać F_VMD (biegun ρ) jest najprostsza.** Zmiana F_V na postać
  Gounarisa–Sakuraia zmienia wyniki o najwyżej około 1 %.
- **Znakowanie fotonem.** Jego naśladowanie w symulacji jest NIESPRAWDZONE.
  Wydajność η znakowanych w danych pochodzi z badania bez cięcia vz.
- **Małe wzorce:** konwersja 307 wierszy, garb 2836 wierszy (rozdział 6.1).
- **Błędy danych.** Błędy danych w przedziałach OA i M(e⁺e⁻) są traktowane jako
  nieskorelowane. Dopasowania wzorców w M4 są słabe (χ²/ndf 2,4–3,6), więc ich
  czynniki porównują kształty, a nie mierzą liczb zdarzeń.
- **Skala masy.** Przeniesienie skali masy wyznaczonej z K⁰_S na M4 jest
  NIESPRAWDZONE.

## 8. Zalecenia dla producentów symulacji

1. **Zdefiniować kanał z modelem pary leptonów.** W `database_decays_eta.dat`
   zapisać parę e⁺e⁻ jako foton wirtualny, tak jak ta sama baza robi to dla
   rozpadu Dalitza π⁰ (`pi0:dilepton[e+,e-],g`). Na przykład:
   `channel=eta:dilepton[e+,e-],pi+,pi-`. Odpowiednio zmienić definicję rozpadu
   w `pluto_decay.h`: zamiast ogólnego „eta --> e+ + e- + pi+ pi-” użyć postaci
   z dileptonem. Czy PLUTO ma wtedy model η → π⁺π⁻γ* (Kroll–Wada, fala P pionów,
   czynniki postaci), jest NIESPRAWDZONE i wymaga testu z punktu 3.
2. **Nie włączać `eta_physics`.** Dodatek jest wyłączony celowo, bo narzuca rozkłady
   z eksperymentu DISTO. Poprawką jest wyłącznie definicja kanału z dileptonem
   z punktu 1.
3. **Przed produkcją sprawdzić generator.** Wystarczy około 10⁵ rozpadów:
   - **M(e⁺e⁻):** udział poniżej 1,52 MeV/c² około 0,07, powyżej 20 MeV/c² około
     0,31; w obecnym generatorze 0,32 i 0,047;
   - **kąt e⁺:** rozkład 1 + cos²θ + (4m_e²/q²) sin²θ, a nie płaski;
   - **M(π⁺π⁻):** z czynnikiem postaci pionu;
   - **zapis:** zapytać, z jaką dokładnością zapisywane są pędy i energie w plikach
     `.evt`. W wydobytej próbce 4,9 % mas pary leży poniżej 2m_e.
4. **Wygenerować zestaw η ponownie** (konwersja, HGeant, DST, PAT). Obok katalogu
   `1_old` są katalogi `smash_eta/1` i `smash_eta/2`. Według producentów
   w listopadzie 2026 r. te same zdarzenia przekonwertowano ponownie, z nowym
   układem kodów generatora. Przed produkcją DST trzeba sprawdzić, jaką definicję
   kanału mają nowe konwersje (NIESPRAWDZONE).
5. **Rozpady η w HGeant.** HGeant wykonuje η → π⁺π⁻π⁰ i η → π⁺π⁻γ według
   przestrzeni fazowej, a η → e⁺e⁻γ z prawdopodobieństwem 0,060 zamiast 0,0069
   (PDG). Ustawia to karta `ETA 37.30 30.20 21.00 4.75 6.0 7.1e-2` w pliku wejściowym
   HGeant (6,0 % zamiast domyślnych w HGeant 0,60 %; producenci zmierzyli 0,064).
   Trzeba ustalić, czy to zwiększenie było celowe. Dla tej analizy lepiej byłoby rozpadać η w PLUTO z modelem wykresu
   Dalitza, a przynajmniej opisać to w dokumentacji produkcji.

Do czasu nowej produkcji obowiązują wagi z rozdziału 5.

## 9. Pliki i makra

Repozytorium `/home/przygoda/HADES/FAT/FAT`, katalog analizy symulacji
`plots/eta_sim/` (indeks: `plots/eta_sim/README.md`, rozdział „Physics
reweighting”).

| plik | zawartość |
|---|---|
| `plots/eta_sim/sim_eta_physics.h` | wzory i źródła, gęstość generatora, wagi (`simphys::RowFactor`), warianty |
| `plots/eta_sim/SimEngine.h` | silnik widm ważonych; opcja `SimOpt::phys` (domyślnie 0 = wyłączona) |
| `plots/eta_sim/sim_rw_build.C` | budowa i sprawdzenie wag (generator, model, sprawdzenie zgodności, kąty, klasy, przekroje czynne) |
| `plots/eta_sim/sim_rw_compare.C` | porównanie z danymi przed i po wagach (pik, garb, selekcje wzbogacone w η, warianty) |
| `plots/eta_sim/sim_rw_dalitz_acc.C` | akceptancja z par z rozpadów Dalitza (krok 2) |
| `plots/eta_sim/sim_eta_peak.C` | prawdziwa M(e⁺e⁻) zaakceptowanego rozpadu bezpośredniego (krok 1) |
| `src/sim_truth.h` | klasy, przynależność, wagi wierszy |

Polecenia (z katalogu repozytorium; katalogi kompilacji poza repozytorium):

```bash
root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_rw_build", kTRUE)' 'plots/eta_sim/sim_rw_build.C+(7)'
root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_rw_compare", kTRUE)' 'plots/eta_sim/sim_rw_compare.C+(15)'
root -l -b -q -e 'gSystem->SetBuildDir("/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/scan/aclic_dalitz_acc", kTRUE)' 'plots/eta_sim/sim_rw_dalitz_acc.C+'
```

Tabele (`plots/eta_sim/tables/`):

- `eta_peak_mee_true.txt`: krok 1;
- `rw_dalitz_acceptance.txt`: krok 2;
- `rw_weights_generator.txt`: gęstość generatora, kąty, dokładność zapisu;
- `rw_weights_model.txt`: wzory, sprawdzenie z PDG, udziały M(e⁺e⁻), średnie wag;
- `rw_weights_closure.txt`: sprawdzenie zgodności na próbce z generatora;
- `rw_weights_bump.txt`: parametry KLOE-2, ⟨f⟩_PS, wB, wC, p3;
- `rw_weights_planes.txt`: człon płaszczyzn (V4/V3);
- `rw_weights_classes.txt`: czynniki wag w klasach dla wariantów;
- `rw_weights_xsec.txt`: przekroje czynne przed i po;
- `rw_weights_peak.txt`: przedziały piku i wydajności cięć dla wariantów;
- `rw_compare_peak.txt`: pik przed i po; rozdział 0 to udziały niezależne od skali;
- `rw_compare_bump.txt`: garb przed i po;
- `rw_compare_enriched.txt`: selekcje wzbogacone w η, dopasowania wzorców,
  metody A–E;
- `rw_compare_variants.txt`: warianty;
- `rw_compare_checks.txt`: 103 sprawdzenia;
- `eta_summary_peak_methods.txt`: metody A–E przed poprawkami.

Duże pliki (`/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/`):

- `generator/eta_pipiee_gen_kin_1000files.txt.gz`: 1 181 801 rozpadów z generatora;
  w wierszu M(e⁺e⁻), M(π⁺π⁻), cos θ e⁺, cos θ π⁺, M(η);
- `scan/rw_generator.root`: gęstość generatora G i stałe;
- `scan/rw_angles_eta.root`: prawdziwe kąty i cos 2φ wierszy klasy 5;
- `scan/sim_compact.root` (std, lep) i `scan/sim_compact_eta.root` (zestaw η):
  pliki pośrednie z prawdziwymi masami;
- `output_{std,lep,eta}_*.root`: wyniki FAT z `scan_nt` z prawdziwymi pędami;
- dzienniki: `scan/rw_v2_build_7.log`, `scan/rw_v2_compare_15_pass1.log`.

Poza tą maszyną:

- konwerter:
  `/lustre/hades/user/hadesdst/sub/feb22/sim/p4500p/convert_smash/smashtrunk/`
  (`database_decays_eta.dat`, `pluto_decay.h`);
- pliki po konwersji: `/lustre/hades/dstsim/feb22/p4500p/smash/convert/smash_eta/1_old/`;
- źródło HGeant: `/home/przygoda/hydra2/hydra2-git/hgeant2-7.4/geant/gdecay.F`.
