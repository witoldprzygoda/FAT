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
