# Analiza π⁺π⁻e⁺e⁻ w p+p przy 4.53 GeV: szybki start

Gałąź `pp45_pippimepem` repozytorium FAT. Wejściem są ntuple PAT z produkcji
„betanew” (29 dni, pliki `lep_034.root` … `lep_067.root`, PID na β_new).
Energia kinetyczna wiązki to 4530 MeV.

## 1. Pobranie kodu

```bash
git clone -b pp45_pippimepem git@github.com:witoldprzygoda/FAT.git FAT
cd FAT
```

Wymagania: ROOT 6 (`root-config` i `hadd` w `PATH`), g++ z obsługą C++17, python3.

## 2. Co trzeba ustawić

**Lista plików wejściowych.** Utwórz listę z własną lokalizacją plików PAT,
jedna ścieżka w wierszu:

```bash
ls /twoja/sciezka/gen4/lep_0*.root > /twoja/sciezka/gen4/lep_betanew.list
wc -l /twoja/sciezka/gen4/lep_betanew.list      # powinno być 29
```

**Konfiguracje.** W trzech plikach `config.json` (π⁺π⁻e⁺e⁻), `config_epep.json`
(π⁺π⁻e⁺e⁺) i `config_emem.json` (π⁺π⁻e⁻e⁻) ustaw listę wejściową (`source`)
i katalog pliku wyjściowego (`filename`); nazwę pliku zostaw bez zmian:

```json
"source": "/twoja/sciezka/gen4/lep_betanew.list",
"filename": "/twoj/katalog/wynikow/output_pippimepem_exp.root",
```

W repozytorium wyniki są skierowane do katalogu
`/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/`. Katalog wyników
musi istnieć przed uruchomieniem (`mkdir -p`). Pozostałych pól nie ruszaj: nazwy
drzew (`PipPimEpEm_ID`, `PipPimEpEp_ID`, `PipPimEmEm_ID`) i energia wiązki
(`kinetic_energy`: 4530 MeV) są już ustawione.

## 3. Kompilacja

```bash
make
```

Powstaje program `./ana`. Kompilację powtórz po każdej zmianie gałęzi lub kodu:
`run_parallel.sh` uruchamia ten `./ana`, który akurat leży w katalogu.

## 4. Uruchomienie

Każdy kanał liczy się osobno. `run_parallel.sh` dzieli listę na części, uruchamia
kilka kopii `./ana` naraz i na końcu scala wyniki przez `hadd`:

```bash
export FAT_TMPDIR=/hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/tmp
mkdir -p $FAT_TMPDIR
./run_parallel.sh config.json      29 15
./run_parallel.sh config_epep.json 29 15
./run_parallel.sh config_emem.json 29 15
```

Argumenty to: liczba części (najwyżej tyle, ile plików na liście) i liczba zadań
uruchamianych naraz, tutaj najwyżej 15. Dwa dalsze są nieobowiązkowe: liczba
ponowień nieudanej części (domyślnie 2; brak gałęzi, wyjątek i przerwanie nie
są ponawiane) i odstęp między startami w sekundach (domyślnie 3). Na maszynie
z mniejszą liczbą rdzeni zmniejsz drugi argument.

- Postęp śledzisz w drugim terminalu: `./monitor_parallel.sh` (z tą samą
  zmienną `FAT_TMPDIR`).
- Pliki pośrednie trafiają do katalogu `FAT_TMPDIR`. Ustawiaj tę zmienną zawsze:
  bez niej skrypt pisze do `/tmp`, a tam brakuje miejsca.
- Część jest nieudana, gdy `./ana` kończy się kodem różnym od zera: w drzewie
  brakuje gałęzi, brakuje pliku z listy, program zgłosił wyjątek w pętli po
  wierszach, wystąpił błąd odczytu lub zapisu (np. pełny dysk) albo zadanie
  przerwano (Ctrl+C). Skrypt wypisuje wtedy numer części, nie scala wyników
  i zostawia katalog roboczy z plikami `cfg_NNN.log`; jego kod wyjścia jest
  różny od zera. Nowy plik wynikowy nie powstaje. Jeśli pod tą nazwą leży plik
  z wcześniejszego uruchomienia, pozostaje bez zmian, a skrypt ostrzega o tym
  i podaje jego datę.
- Wynik odczytaj z komunikatów `run_parallel.sh`, a nie z monitora:
  `monitor_parallel.sh` pokazuje nieudaną część jako trwającą.
- Plik `lep_034.root` w kanale π⁺π⁻e⁺e⁻ liczy się 70 s w jednym zadaniu,
  w kanałach jednoimiennych kilka sekund.

Po zakończeniu w katalogu wyników masz `output_pippimepem_exp.root`,
`output_pippimepep_exp.root` i `output_pippimemem_exp.root`.

## 5. Rysowanie

Makra uruchamiasz z katalogu głównego repozytorium. Wykresy trafiają do
`plots/output/` (PDF i PNG), a całki są wypisywane na konsoli. Makra otwierają
trzy pliki `output_pippim*_exp.root` z katalogu bieżącego, a nie z katalogu
wyników. Przed rysowaniem utwórz więc w katalogu repozytorium dowiązania
symboliczne do nowych plików:

```bash
ln -s /hdd3/przygoda/HADES/FAT_out/exp/pippimepem_isBestAll_T4530/output_pippim*_exp.root .
```

Jeśli w katalogu leżą już starsze pliki o tych nazwach, `ln` zgłosi błąd
i niczego nie nadpisze. Przenieś je wtedy w inne miejsce i powtórz polecenie;
inaczej makra narysują stare dane.

```bash
root -l -b -q plots/mass_spectra_rec_exp.C
root -l -b -q plots/pippimepem_spectra_rec_exp.C
root -l -b -q plots/pippimepem_spectra_cor_exp.C
root -l -b -q plots/pippimepemg_spectra_cor_exp.C
root -l -b -q plots/f1_spectra_cor_exp.C
```

Każdy wykres pokazuje trzy rozkłady: **czarny** to wszystkie kombinacje π⁺π⁻e⁺e⁻,
**czerwony** to tło kombinatoryczne CB = 2√(N₊₊N₋₋) z kanałów jednoimiennych,
**niebieski** to sygnał = wszystkie − CB. Przyrostek `_rec` oznacza pędy
zrekonstruowane, `_cor` pędy z poprawką na straty energii.

| Makro | Co rysuje | Na co patrzeć |
|---|---|---|
| `mass_spectra_rec_exp.C` | M(e⁺e⁻) bez cięcia OA i w oknie 2° ≤ OA ≤ 6° | obszar Dalitza π⁰ (M < 0.14 GeV/c²) i ilość tła CB w oknie OA |
| `pippimepem_spectra_rec_exp.C`, `_cor_exp.C` | M(e⁺e⁻), M(π⁺π⁻), M(π⁺π⁻e⁺e⁻), MM(π⁺π⁻), MM(π⁺π⁻e⁺e⁻); rozkłady kątów używanych w selekcji; M(π⁺π⁻e⁺e⁻) po selekcji i w plasterkach MM | wąski pik η w M(π⁺π⁻e⁺e⁻) w okolicy 0.545 GeV/c²; czy progi selekcji pasują do rozkładów kątów |
| `pippimepemg_spectra_cor_exp.C` | M(π⁺π⁻e⁺e⁻γ) i MM(π⁺π⁻e⁺e⁻γ) z π⁰→e⁺e⁻γ (foton z ECAL, wąskie okno π⁰) | pik ω(782) przy ok. 0.78 GeV/c² i η→π⁺π⁻π⁰ przy ok. 0.55 GeV/c² |
| `f1_spectra_cor_exp.C` | M(π⁺π⁻η) dla η→e⁺e⁻γ i η→γγ | okolica 1.285 GeV/c² (f1); w poprzedniej analizie sygnału nie było |

Uwaga: pik π⁰ w M(γγ) z ECAL leży ok. 15–20% za wysoko, bo energia fotonów nie
jest jeszcze skalibrowana. Okna masy π⁰ mogą więc obejmować tylko część sygnału.

## 6. Własne rozkłady

Główne ntuple to `pippimepem_nt` (pędy zrekonstruowane) i `pippimepem_nt_cor`
(z poprawką strat energii). Mają te same nazwy pól. Każdy wiersz to jedna
zapisana w PAT kombinacja π⁺π⁻e⁺e⁻, która przeszła cięcia: wierzchołek
(`eVertReco_z > −500`) i trigger PT3. Na START nie ma cięcia; `start_iteration`
jest tylko zapisane.

Używane są wszystkie kombinacje (cięcie `isBestAll`, czyli `isBest >= 0`), a nie
tylko najlepsza. Jedno zdarzenie daje więc czasem kilka wierszy, np. tę samą parę
e⁺e⁻ z innym pionem, i każdy wpis w histogramie to kombinacja, a nie zdarzenie.
Wiersze jednego zdarzenia leżą obok siebie i dokładnie jeden z nich ma
`isBest == 1`. Pola `isBest`, `runnumber` i `event` są w ntuplach, ale para
(`runnumber`, `event`) nie wyznacza zdarzenia jednoznacznie.

Ta sama para leptonów może wystąpić w kilku wierszach zdarzenia, więc wielkości
zależne tylko od niej (`m_ee`, `oa`) wchodzą do rozkładów kilka razy: w pliku
`lep_034.root` wierszy e⁺e⁻ jest 83424, a zdarzeń 70394. Do widm samej pary
leptonów, także w kanałach jednoimiennych, użyj warunku `isBest == 1` albo licz
każdą parę raz na zdarzenie.

Ntupla `epemg_nt` ma jeden wiersz na kombinację i na foton z ECAL, który
przeszedł cięcia jakości, także dla kombinacji z `isBest == 0`. Pól `isBest`,
`runnumber` i `event` w niej nie ma. Jej wiersze idą w kolejności `scan_nt`
(rozdział 7): wierszowi `scan_nt` z `n_gam` = k odpowiada k kolejnych wierszy
`epemg_nt`.

Flagi do własnych selekcji:

| Flaga | Znaczenie |
|---|---|
| `oa_pass` | 2° ≤ OA(e⁺e⁻) ≤ 6° (sam kąt jest w polu `oa`) |
| `sel_pass` | wszystkie trzy naraz: OA_LAB((π⁺π⁻),(e⁺e⁻)) < 50°, M(π⁺π⁻) < 0.42 GeV/c², OA w układzie η > 140° |
| `pippimepemg_pass_narrow` | jeden foton z ECAL i M(e⁺e⁻γ) w [0.125, 0.145] GeV/c² |
| `eta_dalitz_pass` | jeden foton z ECAL i M(e⁺e⁻γ) w [0.5, 0.6] GeV/c² |
| `eta_gg_pass` | dwa fotony z ECAL i M(γγ) w [0.5, 0.6] GeV/c² |

Tło CB dla dowolnego rozkładu liczysz bin po binie z tych samych ntupli w
plikach `output_pippimepep_exp.root` i `output_pippimemem_exp.root`:
CB = 2√(N₊₊·N₋₋). Gotowa funkcja to `PlotUtils::drawSignal` w `plots/PlotUtils.h`,
którą można wywołać na dowolnej zmiennej i z dowolnym warunkiem.

## 7. Ntupla `scan_nt` do badania cięć

`scan_nt` ma te same wiersze w tej samej kolejności co `pippimepem_nt`: jeden
wiersz na kombinację po cięciach `isBestAll`, wierzchołka i triggera. Innej
wstępnej selekcji nie ma. Każde cięcie jest zapisane dwa razy: jako wartość
zmiennej i jako flaga 1/0 policzona dla każdego wiersza, niezależnie od
pozostałych cięć. Nazwa bez przyrostka oznacza pędy zrekonstruowane, przyrostek
`_cor` pędy z poprawką.

| Grupa pól | Zawartość |
|---|---|
| klucze i zdarzenie | `file_id` (numer z nazwy pliku PAT, `lep_034.root` → 34), `pat_entry` (numer wiersza w drzewie tego pliku PAT), `runnumber`, `event`, `isBest`, `trigbit`, `hyp_chi2` (χ² hipotezy; wiersz z najmniejszą wartością w zdarzeniu ma `isBest == 1`), `eVertReco_x/y/z/chi2`, `neutr_mult`, `fwdet_mult`, `evtStart`, `start_iteration` |
| ślady `ep_`, `em_`, `pip_`, `pim_` | `p_rec`, `p_cor`, `theta`, `phi`, `beta_new`, `system`, `sector`, `r`, `z`, `rkchi2`, `isUsedVertex`, `oa_lept`, `oa_hadr`, `tof_new`, `tof_mom` |
| tylko leptony | `theta_rich`, `phi_rich`, `rich_ind`, `richmatchquality`, `rich_padnum`, `ringchi2`, `kIsLepton` |
| tylko piony | `isring`, `dedx_mdc`, `mass2` = p_rec²(1/β_new² − 1) w (MeV/c²)² |
| zmienne cięć | `oa`; z przyrostkiem `_cor` i bez: `m_pippimepem`, `m_ee`, `m_pippim`, `oa_pippim_epem_lab`, `oa_pippim_epem_eta_rest`, `mm_pippimepem`, `mm_pippim`, `mm_epem`, `pt4` (pęd poprzeczny układu π⁺π⁻e⁺e⁻, MeV/c), `y4_cms` (pospieszność tego układu w układzie środka masy wiązki i tarczy) |
| flagi | `pass_oa_min` (OA(e⁺e⁻) > 2°), `pass_oa_max` (OA(e⁺e⁻) < 6°), `oa_pass`; z przyrostkiem `_cor` i bez: `pass_oa_lab` (< 50°), `pass_m_pippim` (< 0.42 GeV/c²), `pass_oa_eta_rest` (> 140°), `sel_pass` (wszystkie trzy naraz); `cor_ok` (1, gdy wszystkie cztery poprawione pędy są liczbami, 0, gdy któryś jest NaN) |
| ECAL | `n_gam` (liczba klastrów neutralnych, które przeszły cięcia jakości) i trzy z nich o największej energii: `g1_e`, `g1_theta`, `g1_phi`, … `g3_phi` |
| detektor przedni | dwaj pierwsi kandydaci: `fw1_p`, `fw1_theta`, `fw1_phi`, `fw2_…` |

Brak fotonu lub kandydata w detektorze przednim jest oznaczony wartością −1;
w pozostałych polach −1 nie ma takiego znaczenia. Ujemna masa brakująca oznacza
ujemny kwadrat masy (zapis −√(−M²)). Pole `mass2` ma wartość −10⁶, gdy β_new ≤ 0.

`neutr_mult` i `fwdet_mult` to krotności z PAT; brak klastra lub kandydata
(w PAT −100) jest zapisany jako 0. PAT zapisuje najwyżej 5 klastrów neutralnych
i 3 kandydatów w detektorze przednim, a krotność bywa większa (w `lep_034.root`
do 8 i do 5). `n_gam`, `g1`–`g3`, `fw1` i `fw2` dotyczą tylko zapisanych obiektów.

Kąt `phi` śladów jest w zakresie od −180° do 180°, a `phi_rich`, `g*_phi` i
`fw*_phi` od 0° do 360°, jak w PAT; przed odjęciem dwóch kątów sprowadź je do
wspólnego zakresu. Wiersz PAT odpowiadający wierszowi `scan_nt` wskazuje para
(`file_id`, `pat_entry`), niezależnie od podziału listy na części.

W części wierszy PAT nie podaje poprawionego pędu π⁻ (`pim_p_cor` to NaN; są to
π⁻ pod kątem θ < 16°). Wielkości `_cor` z udziałem π⁻ są wtedy NaN, flagi `_cor`
mają wartość 0 i `cor_ok` = 0; pęd zrekonstruowany nie jest podstawiany.
W `lep_034.root` dotyczy to 0.7% wszystkich wierszy e⁺e⁻, ale 1.9% wierszy
z `sel_pass == 1` i 2.9% takich wierszy w oknie η (0.53 < M(π⁺π⁻e⁺e⁻) < 0.56
GeV/c²). Licząc cokolwiek na polach `_cor`, dodaj warunek `cor_ok == 1` i podaj,
ile wierszy odpadło. Warunek zaprzeczony, np. `!(m_pippim_cor >= 0.42)`,
przepuszcza wiersze z NaN.

## 8. Symulacja

Ta sama analiza działa na symulacji SMASH: na zbiorze standardowym i na zbiorze
leptonowym (zbiór η jest przewidziany w kodzie, ale jego plików jeszcze nie ma).
Wszystko, co dotyczy symulacji, włączają klucze konfiguracji, których wartości
domyślne odpowiadają danym, więc wyniki dla danych się nie zmieniają.

| Konfiguracja | Zbiór (lista) | Kanał |
|---|---|---|
| `config_sim_std.json`, `config_sim_std_epep.json`, `config_sim_std_emem.json` | standardowy (`sim_std.list`) | π⁺π⁻e⁺e⁻, π⁺π⁻e⁺e⁺, π⁺π⁻e⁻e⁻ |
| `config_sim_lep.json`, `config_sim_lep_epep.json`, `config_sim_lep_emem.json` | leptonowy (`sim_lep.list`) | jak wyżej |

Każda lista zawiera jedną ścieżkę pliku. Wyniki trafiają do
`/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/output_<zbiór>_<kanał>.root`.
Względem konfiguracji danych zmieniają się trzy pola:

- `"kinematics": {"simulated": true}` to główny przełącznik: program czyta
  prawdę z GEANT o czterech śladach i zapisuje ją w `scan_nt`;
- `"simulation": {"set": 1, "pdg_br": true, "br_eta_pipiee": 2.68e-4, "br_x17": 0.0}`
  podaje numer zbioru (1 standardowy, 2 leptonowy, 3 η) i stosunki
  rozgałęzień użyte w wagach;
- `"trigger": {"bit_test": true}`: PT3 oznacza tu ustawiony bit 8192, czyli
  `(trigbit & 8192) != 0`. W symulacji wiersze PT3 mają `trigbit` = 12288;
  w danych zostaje warunek `trigbit == 8192`.

Drzewa symulacji nie mają pola `start_iteration`; w ntuplach ma ono wartość −1.
Pole `pat_entry` jest typu float i dokładnie zapisuje numery wierszy tylko do
2²⁴ = 16 777 216. Dlatego w symulacji program odmawia pracy (kod 1, bez pliku
wynikowego), gdy plik wejściowy ma więcej wierszy; taki plik trzeba podzielić.

**Uruchomienie.** Gdy lista zawiera jeden plik, a części ma być więcej niż jedna,
`run_parallel.sh` dzieli wiersze drzewa tego pliku na równe przedziały kolejnych
numerów wierszy:

```bash
export FAT_TMPDIR=/hdd3/przygoda/HADES/FAT_out/sim/pippimepem_betanew/tmp
./run_parallel.sh config_sim_std.json 4 4
./run_parallel.sh config_sim_lep.json 15 15
```

Skrypt wypisuje przedział każdej części. Przed scaleniem sprawdza w dziennikach
zadań, czy każda część przetworzyła dokładnie swój przedział i czy przedziały
pokrywają całość bez przerw i bez nakładania się. Jeśli nie, nie scala wyników
i kończy się kodem 5. Pola `start_event` i `max_events` zawężają dzielony
przedział (przydaje się do prób). Z pliku leptonowego czytany jest cykl `;3`
drzewa `PipPimEpEm_ID` (7 622 819 wierszy), bo łańcuch plików bierze cykl
najwyższy; cykl `;2` to tylko jego pierwsze 55 454 wiersze.

Podziału na przedziały wierszy nie ma (plik idzie w jednej części, skrypt
wypisuje ostrzeżenie), gdy:

- `trigger.bias_correction` jest włączone: program liczy wagę PT3 pliku
  z wierszy, które czyta jedno zadanie, więc podział zmieniłby `trigger_corr`;
- wiersz listy nie jest zwykłym plikiem (wzorzec z `*`, adres sieciowy).

Inne zasady skryptu:

- Liczba części może wynosić od 1 do 999 (numery trzycyfrowe utrzymują
  kolejność scalania). Bez trzeciego argumentu naraz działa najwyżej 15 zadań.
- Przed scaleniem skrypt sprawdza kolumny ntupli w częściach. Część, w której
  ntupla nie ma ani jednego wiersza, zapisuje ją tylko z kolumną `trigger_corr`.
  Gdyby taka część szła pierwsza, `hadd` zgubiłby pozostałe kolumny bez żadnego
  błędu. Dlatego skrypt zapisuje takie puste ntuple od nowa z kolumnami części,
  które mają wiersze. Po scaleniu sprawdza, czy każda ntupla ma te kolumny
  i sumę wierszy części. Jeśli nie, kończy się kodem 3 i nie zapisuje wyniku.
- Przerwanie skryptu (Ctrl+C, `kill`, zamknięta sesja ssh) zatrzymuje
  wszystkie działające części; nic nie jest scalane, a katalog roboczy
  zostaje (kody 130, 143, 129). Długie przebiegi najlepiej uruchamiać przez
  `nohup` z wyjściem do pliku dziennika.

**Dodatkowe pola `scan_nt`.** Do 145 pól danych dochodzi 88 pól:

| Pola | Zawartość |
|---|---|
| `ep_`, `em_`, `pip_`, `pim_` + `sim_id`, `sim_parentid`, `sim_grandparentid`, `sim_processid`, `sim_mediumid`, `sim_p`, `sim_px`, `sim_py`, `sim_pz`, `sim_vertexx`, `sim_vertexy`, `sim_vertexz`, `sim_geninfo`, `sim_geninfo1`, `sim_geninfo2`, `sim_genweight`, `sim_iscommon`, `sim_isghost`; leptony także `rich_sim_geanttrack` | prawda z GEANT o śladzie, przepisana z PAT (w kanałach jednoimiennych `ep_` i `em_` to pierwszy i drugi lepton) |
| `sim_set`, `sim_rules_ver` | numer zbioru i wersja reguł (2) |
| `sim_cls`, `sim_lcls`, `sim_pcls` | klasa czwórki, pary leptonów (L0–L11) i pary pionów (P0–P9) |
| `sim_w`, `sim_w_db` | waga wiersza ze stosunkami rozgałęzień PDG i z bazy konwertera |
| `sim_own_A`, `sim_own_B` | 1, gdy wiersz należy do tego zbioru w podziale A (podstawowym) lub B (wariant) |
| `sim_flags` | bity ostrzeżeń (tabela niżej) |
| `sim_tmpl` | wzorzec: 1 sygnał, 2 tło fizyczne, 3 fałszywy lepton, 4 leptony z różnych źródeł (opis niżej) |
| `m4_true`, `mee_true`, `mpipi_true` | prawdziwe M(π⁺π⁻e⁺e⁻), M(e⁺e⁻), M(π⁺π⁻) w GeV/c² z pędów GEANT i mas nominalnych; −1, gdy któryś ślad jest inną cząstką |

Reguły klas, przynależności i wag są w `src/sim_truth.h` (wersja reguł 2 z 13
listopada 2026 r.; wersja 1 była liczbowo taka sama jak wzorcowy `simrules.h`,
a wersja 2 różni się od niej tylko stosunkami rozgałęzień rozpadów π⁰ i η
wykonywanych przez HGeant). Te stosunki rozgałęzień pochodzą teraz wprost z karty
`ETA 37.30 30.20 21.00 4.75 6.0 7.1e-2` w plikach wejściowych HGeant tej produkcji
(dla π⁰ z wartości domyślnych HGeant), a nie z pomiaru producentów na 10⁶
zdarzeniach (raport `docs/latex/report.tex`, rozdział „The HGeant decay table of
the eta”). Plik nie zależy od ROOT, więc
można go dołączyć w makrze i przeliczyć wagi z innymi stosunkami rozgałęzień.

Klasy czwórki (`sim_cls`): 0 błędnie zidentyfikowana cząstka, 1 pion rozpadł
się na mion, 2 garb (η → π⁺π⁻π⁰, π⁰ → e⁺e⁻γ), 3 garb z konwersją fotonu
z π⁰, 4 pik z konwersją fotonu (η → π⁺π⁻γ), 5 pik η → π⁺π⁻e⁺e⁻ (zbiór η),
6 X17 (zbiór η), 7 piony z parą e⁺e⁻ z jednego źródła, 8 leptony z różnych
źródeł.

**Wagi i przynależność.** Do widm bierz wiersze z `sim_own_A == 1` (wariant:
`sim_own_B == 1`) z wagą `sim_w`; wiersze z zerem dostarcza drugi zbiór.
Tych samych rodzajów wierszy nigdy nie dodawaj z obu zbiorów. Oczekiwana
liczba wierszy w danych to (L_dane / 221.05 nb⁻¹) · Σ `sim_w`, przy czym
L_dane jest parametrem. Wariancja w przedziale histogramu to suma po
zdarzeniach kwadratów sum wag ich wierszy; zdarzenie to ciąg sąsiednich
wierszy o tych samych `event` i `eVertReco_z` (w symulacji `runnumber` = 0,
`file_id` = −1). Histogramy w plikach symulacji zawierają tylko wiersze
z `sim_own_A == 1` i nie są ważone (`trigger_corr` = 1), więc nadają się
jedynie do kontroli. Widma licz z `scan_nt`, w której są wszystkie wiersze.
Cięcie `eVertReco_z > −160` nakładasz sam, jak w danych.

`sim_w_db` to wariant z wagami z bazy konwertera. Nie zawiera poprawek
rozpadów naturalnych: zostają w nim stosunki rozgałęzień HGeant. Na przykład
para z rozpadu Dalitza η w zbiorze standardowym ma `sim_w_db` = 1, a nie
0.114 (w e⁺e⁻, wiersze należące do zbioru standardowego, vz > −160:
Σ `sim_w` = 371, Σ `sim_w_db` = 3250). Wariant ma więc sens tylko dla
wierszy z wymuszonymi rozpadami, które należą do zbioru leptonowego.

**Sygnał bez tła kombinatorycznego.** W symulacji sygnałem kanału są wiersze,
których cząstki według kodów generatora pochodzą dokładnie z tego kanału
produkcji i rozpadu (`sim_cls` od 2 do 6). Niczego się nie odejmuje, a drzewa
jednoimienne nie służą za tło, tylko do badań pobocznych (isBest,
kombinatoryka); czynnik k pochodzi wyłącznie z danych.

Podział na wzorce podaje pole `sim_tmpl`:

| `sim_tmpl` | Wzorzec | Warunek |
|---|---|---|
| 1 | sygnał | `sim_cls` od 2 do 6, także z leptonem-duchem |
| 2 | tło fizyczne | prawdziwa para e⁺e⁻ z jednego źródła (L1–L6; w zbiorze η także L10, L11 poza sygnałem), dowolne piony (również błędnie zidentyfikowane i rozpadłe na mion), lepton nie jest fałszywy |
| 3 | fałszywy lepton | lepton błędnie zidentyfikowany (L0), ślad-duch na miejscu leptonu (`sim_isghost` = 1) albo ta sama cząstka na obu miejscach leptonów |
| 4 | leptony z różnych źródeł | L7–L9; ani sygnał, ani tło (w symulacji nie ma tła kombinatorycznego) |

`sim_cls` 0, 1 i 7 to nie są te wzorce: klasy 0 i 1 zawierają też wiersze
z prawdziwą parą e⁺e⁻ i fałszywym pionem, a klasa 7 także wiersze
z leptonem-duchem. Wielkości w zbiorze standardowym (e⁺e⁻, wiersze należące
do zbioru w podziale A, vz > −160, 115 414 wierszy):

| `sim_tmpl` | Wiersze | Σ `sim_w` |
|---|---|---|
| 1 sygnał | 3 143, w tym 395 z leptonem-duchem | 3309, w tym 413 |
| 2 tło fizyczne | 86 834, w tym 13 572 z `sim_cls` 0 lub 1 | 83 999 |
| 3 fałszywy lepton | 22 370: 1 401 L0 i 20 969 z leptonem-duchem | 22 063 |
| 4 leptony z różnych źródeł | 3 067 | 2643 |

Wiersze sygnału z leptonem-duchem (12.6%) zostają w sygnale; bit 10 pola
`sim_flags` pozwala je wyłączyć. W drzewach jednoimiennych wszystkie wiersze
z parą z jednego źródła (L1–L6) mają tę samą cząstkę na obu miejscach
leptonów, więc trafiają do wzorca 3.

| Bit `sim_flags` | Znaczenie |
|---|---|
| 0 (1) | zapisana waga generatora produktu wymuszonego rozpadu różni się od tabeli |
| 1 (2) | zarezerwowany, nigdy nie jest ustawiany (rozpad wymuszony bez kodu w tabeli wag nie może wystąpić) |
| 2 (4) | produkt drugiego rozpadu (π⁰ z ω → π⁰e⁺e⁻) razem z leptonem z dowolnego rozpadu ω → π⁰e⁺e⁻ w wierszu (licznika rozpadu się nie porównuje); czynnik łańcucha nie jest doliczony |
| 3 (8) | kod generatora spoza słownika |
| 4 (16) | zaokrąglony kod pasuje do kilku kodów słownika (sprawdzane są wszystkie liczby całkowite, które po zapisie jako float dają tę wartość) |
| 5 (32) | elektron lub pozyton na miejscu pionu |
| 6 (64) | w wierszu jest produkt drugiego rozpadu (525107 / 7051) |
| 7 (128) | więcej niż jeden wymuszony rozpad (waga jest iloczynem) |
| 8 (256) | ślad z wagą > 0 bez wymuszonego kodu tego zbioru |
| 9 (512) | ślad z wymuszonym kodem innego zbioru (zły `simulation.set`?) |
| 10 (1024) | ślad-duch na miejscu leptonu |
| 11 (2048) | ślad-duch na miejscu pionu |
| 12 (4096) | ta sama cząstka GEANT na obu miejscach leptonów (ten sam rodzaj, pęd i punkt powstania); w drzewie e⁺e⁻ jedno z miejsc ma wtedy zły rodzaj cząstki, więc wiersz i tak jest w L0 |
| 13 (8192) | zbiór leptonowy: produkt wymuszonego rozpadu ρ⁰/ω/φ razem ze śladem FC (np. konwersja fotonu z rozpadu Dalitza π⁰ z ω → π⁰e⁺e⁻); według reguł wersji 1 i 2 takiego wiersza nie dostarcza żaden zbiór |
| 14 (16384) | kod, którego reguły wersji 1 i 2 nie obsługują: odmiany φ → K⁰K̄⁰ 551010, 551016, 551610 (brak ich na liście rozpadów hadronowych, więc wiersz zbioru standardowego z leptonem FD jest odrzucany w podziale A, choć nie dostarcza go żaden zbiór) albo produkt η z rozpadu φ → ηγ (550117, bez poprawki η) |

Bity 13 i 14 oznaczają znane luki reguł wersji 1, które zostały też w wersji 2. Ich zmiana zmieniłaby
przynależność i wagi wzorcowych tabel, więc reguły zostały bez zmian,
a wiersze są tylko oznaczone. W zbiorze standardowym (e⁺e⁻) bit 14 ma 27
wierszy, z czego 1 jest odrzucony w podziale A. W zbiorze leptonowym (e⁺e⁻,
wiersze 0–300 000, vz > −160) bit 13 ma 338 wierszy; to około 1% wagi
wierszy z rozpadami ρ⁰/ω/φ.

Rozkład M(e⁺e⁻) rozpadu η → π⁺π⁻e⁺e⁻ w zbiorze η jest niezgodny z QED już w generatorze; obserwację i poprawki w postaci wag opisuje [SYMULACJA_ETA_GENERATOR.md](SYMULACJA_ETA_GENERATOR.md).
