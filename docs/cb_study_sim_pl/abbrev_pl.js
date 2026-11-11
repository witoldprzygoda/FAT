// abbrev_pl.js -- skróty i symbole, rozwijane na slajdzie, na którym pojawiają się po raz pierwszy.
// Klucz: dokładny tytuł slajdu. Tekst trafia do przypisu i do notatek (lib_pl.js). Znaczniki jak w rt(): **pogrubienie**, _{dolny}, ^{górny}.
module.exports = {
  'Tło kombinatoryczne e⁺e⁻ w symulacji SMASH: sygnał, pary jednoimienne i wybór isBest':
    '**HADES**: High Acceptance Di-Electron Spectrometer · **p+p**: zderzenia proton–proton · **T**: energia kinetyczna wiązki · **SMASH**: Simulating Many Accelerated Strongly-interacting Hadrons, generator zderzeń · **gen4**: czwarta generacja rekonstrukcji · **isBest**: znacznik wyboru pary w zdarzeniu',
  'W skrócie: wzór na tło działa dla wszystkich kombinacji, a przy isBest == 1 zawodzi':
    '**N_{++}, N_{−−}**: liczby par e⁺e⁺ i e⁻e⁻ w przedziale masy · **M_{ee}**: masa niezmiennicza pary · **isBest ≥ 0**: wszystkie zapisane kombinacje · **isBest == 1**: jedna wybrana para w zdarzeniu',
  'Próbka: 8,04 mln par e⁺e⁻ z symulacji SMASH i pełna historia każdego śladu':
    '**GEANT**: program symulujący przejście cząstek przez detektor · **PAT**: PostDST Analysis Tool, program zapisujący hipotezy par · **β_{new}**: prędkość wyznaczona z czasu względnego · **PID**: identyfikacja cząstek',
  'Kanały: pary prawdziwe pochodzą z rozpadów Dalitza π⁰ i η oraz z konwersji fotonów':
    '**OA**: kąt otwarcia pary · **γ**: foton · **rozpad Dalitza**: rozpad mezonu na foton i parę e⁺e⁻ · **konwersja**: zamiana fotonu na parę e⁺e⁻ w materiale detektora',
  'Pochodzenie π⁰ w rozpadach Dalitza: około 90 % zapisano jako pierwotne, około 6 % pochodzi z ρ i ω':
    '**kod geninfo**: liczba zapisana w symulacji, opisująca przodka śladu i produkty jego rozpadu; dwukrotnie powtórzony numer cząstki oznacza cząstkę pierwotną (707 to π⁰, 1717 to η)',
  'Identyfikacja leptonów: β_{new} względem q·p, cięcie 5σ':
    '**q**: znak ładunku śladu · **p**: pęd · **σ**: szerokość pasma leptonów w zmiennej β_{new}, wyznaczona osobno dla symulacji',
  'Selekcja: wszystkie użyte cięcia':
    '**eVertReco_z**: współrzędna z zrekonstruowanego wierzchołka zdarzenia · **OA**: kąt otwarcia pary · **RICH**: detektor promieniowania Czerenkowa (ang. Ring Imaging Cherenkov)',
  'Cztery wielkości na widmie masy e⁺e⁻ i sposób ich liczenia':
    '**CB**: tło kombinatoryczne (ang. combinatorial background) · **N_{+−}**: liczba par e⁺e⁻ w przedziale masy · **S**: sygnał · **k**: poprawka akceptancji par jednoimiennych względem różnoimiennych',
  'Co robi isBest: w zdarzeniu rywalizują tylko hipotezy tego samego rodzaju':
    '**hipoteza**: jedna kombinacja dwóch leptonów w zdarzeniu · **χ²**: miara zgodności zmierzonych czasów przelotu z hipotezą pary',
  'Dlaczego isBest nie wskazuje pary prawdziwej: kryterium χ² porównuje tylko czasy przelotu':
    '**RPC, TOF**: detektory czasu przelotu HADES · **χ²**: miara zgodności czasów przelotu z hipotezą pary · **PAT3**: wersja programu PAT użyta w tej produkcji',
  'Jak wyznaczyć poprawkę k: mieszanie zdarzeń w danych':
    '**k_{mix}**: poprawka z mieszania zdarzeń · **PT3**: tryger fizyczny 3 · **n**: liczba leptonów w zdarzeniu · **F_{n}**: wzór 2√(N_{++} N_{−−}) w klasie n',
  'Do czego prowadzi isBest ≥ 0: wzór zgadza się z tłem, a istotność spada tylko o 2 %':
    '**S/√(N + σ²)**: istotność statystyczna; S to sygnał prawdziwy, N liczba par e⁺e⁻, σ² wariancja tła ze wzoru',
  'Wnioski: wszystkie kombinacje do widma inkluzywnego; isBest == 1 traci sygnał i wymaga poprawki tła':
    '**PT3**: tryger fizyczny 3 (ang. physics trigger) · **klasa n**: liczba zidentyfikowanych leptonów w zdarzeniu',
};
