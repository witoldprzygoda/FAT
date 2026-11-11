# Tło kombinatoryczne e⁺e⁻ w symulacji SMASH: prezentacja

Talia po polsku, w stylu opisanym w `/hdd3/przygoda/HADES/files/sim/SLIDES_STYLE.md`.

| Plik | Rola |
|---|---|
| `Tlo_kombinatoryczne_SMASH_PL.pptx` | gotowa prezentacja |
| `export_deck_data.C` | makro ROOT: z `plots/output/cb_study_sim_std_full.root` i plików PID tworzy `data/deck_numbers.json` |
| `export_extra.C` | makro ROOT: dodatkowe histogramy (kontury PID v6, mieszanie zdarzeń, dane, test wyboru isBest) do `data/deck_extra.json`; wejście: `plots/k_mix_pools.C`, `plots/isbest_choice_sim.C`, `plots/isbest_fraction_betanew.C`, `plots/isbest_emulation_betanew.C`, `mixing_betanew.root` |
| `make_assets.py` | rysunki (`figs/`), równania (`eq/`, wymaga pdflatex i mutool), liczby do tabel (`data/deck_tables.json`) |
| `lib_pl.js` | kopia `lib.js` ze stopką, tytułem i typografią polską |
| `abbrev_pl.js` | przypisy ze skrótami; klucz to dokładny tytuł slajdu |
| `build_deck.js` | treść slajdów i notatki |
| `rebuild.sh` | budowanie, walidacja, obrazy slajdów |

Budowanie:

```
TOOLS=<katalog z node_modules/pptxgenjs> SKILL=<katalog skilla pptx> ./rebuild.sh          # sama talia
TOOLS=<katalog z node_modules/pptxgenjs> SKILL=<katalog skilla pptx> ./rebuild.sh assets   # także rysunki i liczby
```

`build_deck.js` po zapisie usuwa powtórzone znaczniki `<a:pPr>`, które pptxgenjs wstawia przed każdym kolejnym
fragmentem akapitu: bez tego akapity z pogrubieniem albo indeksem tracą punktor.

Kontury PID na rysunku pochodzą z pliku v6, a liczby z produkcji z konturami v5 (różnice promilowe, bez przeliczania). Poprawka k pochodzi wyłącznie z mieszania zdarzeń w danych: `plots/k_mix_pools.C` i `make_assets.py`; symulacja służy tylko do testu metody. Test zależności od składu puli: `diag/k_pool_composition_sim.cc`.

Tabele biorą liczby z `data/deck_tables.json`. Liczby wpisane wprost w `build_deck.js` pochodzą z
`plots/output/cb_study_sim_std_full.txt` i `plots/output/cb_study_emu_sim_std_full.txt`.
