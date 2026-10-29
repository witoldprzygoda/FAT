# Analiza π⁺π⁻e⁺e⁻ w p+p przy 4.5 GeV: szybki start

Gałąź `pp45_pippimepem` repozytorium FAT. Wejściem są ntuple PAT z produkcji
„betanew” (29 dni, pliki `lep_034.root` … `lep_067.root`, PID na β_new).

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
(π⁺π⁻e⁺e⁺) i `config_emem.json` (π⁺π⁻e⁻e⁻) zmień tylko pole `source`:

```json
"source": "/twoja/sciezka/gen4/lep_betanew.list",
```

Pozostałych pól nie ruszaj. Nazwy drzew (`PipPimEpEm_ID`, `PipPimEpEp_ID`,
`PipPimEmEm_ID`) i plików wyjściowych są już ustawione. Pliki wyjściowe powstają
w katalogu głównym repozytorium.

## 3. Kompilacja

```bash
make
```

Powstaje program `./ana`.

## 4. Uruchomienie

Każdy kanał liczy się osobno. `run_parallel.sh` dzieli listę na części, uruchamia
kilka kopii `./ana` naraz i na końcu scala wyniki przez `hadd`:

```bash
./run_parallel.sh config.json      29 29 2 4
./run_parallel.sh config_epep.json 29 29 2 4
./run_parallel.sh config_emem.json 29 29 2 4
```

Argumenty to: liczba części (najwyżej tyle, ile plików na liście), liczba zadań
uruchamianych naraz, liczba ponowień nieudanej części i odstęp między startami w
sekundach. Na maszynie z mniejszą liczbą rdzeni zmniejsz drugi argument, np. `29 8 2 4`.

- Postęp śledzisz w drugim terminalu: `./monitor_parallel.sh`.
- Pliki pośrednie trafiają do `/tmp`. Jeśli brakuje tam miejsca, wskaż inny katalog:
  `FAT_TMPDIR=/duzy/dysk ./run_parallel.sh ...`
- Na 32 rdzeniach kanał π⁺π⁻e⁺e⁻ liczy się kilka minut.

Po zakończeniu masz `output_pippimepem_exp.root`, `output_pippimepep_exp.root`
i `output_pippimemem_exp.root`.

## 5. Rysowanie

Makra uruchamiasz z katalogu głównego repozytorium. Wykresy trafiają do
`plots/output/` (PDF i PNG), a całki są wypisywane na konsoli.

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
kombinacja z `isBest==1`, która przeszła cięcia: wierzchołek (`eVertReco_z > −500`)
i trigger PT3. Na START nie ma cięcia; `start_iteration` jest tylko zapisane.

Flagi do własnych selekcji:

| Flaga | Znaczenie |
|---|---|
| `oa_pass` | 2° ≤ OA(e⁺e⁻) ≤ 6° (sam kąt jest w polu `oa`) |
| `sel_pass` | OA_LAB((π⁺π⁻),(e⁺e⁻)) < 50°, M(π⁺π⁻) < 0.42 GeV/c², OA w układzie η > 140° |
| `pippimepemg_pass_narrow` | jeden foton z ECAL i M(e⁺e⁻γ) w [0.125, 0.145] GeV/c² |
| `eta_dalitz_pass` | jeden foton z ECAL i M(e⁺e⁻γ) w [0.5, 0.6] GeV/c² |
| `eta_gg_pass` | dwa fotony z ECAL i M(γγ) w [0.5, 0.6] GeV/c² |

Tło CB dla dowolnego rozkładu liczysz bin po binie z tych samych ntupli w
plikach `output_pippimepep_exp.root` i `output_pippimemem_exp.root`:
CB = 2√(N₊₊·N₋₋). Gotowa funkcja to `PlotUtils::drawSignal` w `plots/PlotUtils.h`,
którą można wywołać na dowolnej zmiennej i z dowolnym warunkiem.
