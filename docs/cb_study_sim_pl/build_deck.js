// build_deck.js -- prezentacja o tle kombinatorycznym e⁺e⁻ (symulacja SMASH). Styl: SLIDES_STYLE.md, biblioteka lib_pl.js.
// Liczby do tabel: data/deck_tables.json (make_assets.py). Uruchamianie: ./rebuild.sh
const L = require('./lib_pl.js');
const { fs, path, pres, C, MX, CW, W, rt, tx, bl, paras, section, cs, fig, cap, eq, result, card, stat, tbl, flushFootnotes, withAbbr, footLater } = L;
const T = JSON.parse(fs.readFileSync(path.join(__dirname, 'data/deck_tables.json')));
const OUT = 'Tlo_kombinatoryczne_SMASH_PL.pptx';
const M = 'GeV/c^{2}';

function titleSlide(sec, title, sub, notes) {
  section(sec);
  const s = pres.addSlide({ masterName: 'TITLE', sectionTitle: sec });
  s.addText(rt(title), { placeholder: 'title' });
  s.addText(rt(sub), { placeholder: 'subtitle' });
  s.addNotes(withAbbr(title, notes));
  footLater(s, title);
}
function divider(sec, title, items, notes) {
  section(sec);
  const s = pres.addSlide({ masterName: 'DIVIDER', sectionTitle: sec });
  s.addText(rt(title), { placeholder: 'title' });
  bl(s, items, { x: MX, y: 2.35, w: CW, h: 3.0, fontSize: 20 });
  s.addNotes(notes);
}
const num = { align: 'right' };
const pct = t => (parseFloat(t.replace(',', '.'))*100).toFixed(1).replace('.', ',') + ' %';
const D4 = T.data[4], num2 = v => (100*(1 - v)).toFixed(0), nv = t => parseFloat(t.replace(/\s/g, ''));
const cbTooBig = ((nv(D4.Cb)/nv(D4.Ce) - 1)*100).toFixed(0);
const P1 = n => `Część 1 z 5 · ${n}`, P2 = n => `Część 2 z 5 · ${n}`, P3 = n => `Część 3 z 5 · ${n}`, P4 = n => `Część 4 z 5 · ${n}`, P5 = n => `Część 5 z 5 · ${n}`;
const A = { fill: C.accent1, color: C.accent1 }, B_ = { fill: C.accent2, color: C.accent2 };   // isBest >= 0 blue, isBest == 1 orange

// ================================================================ wprowadzenie
titleSlide('Wprowadzenie',
  'Tło kombinatoryczne e⁺e⁻ w symulacji SMASH: sygnał, pary jednoimienne i wybór isBest',
  'Reakcja p+p przy T = 4,53 GeV. Symulacja SMASH standard (gen4) z pełną historią śladów.',
  'Studium odpowiada na trzy pytania. Jak dobrze wzór 2√(N++ N−−) odtwarza tło kombinatoryczne policzone dokładnie ze śladów w symulacji? Co dzieje się z tym wzorem, gdy z każdego zdarzenia bierze się jedną parę (isBest == 1)? Do czego prowadzi użycie wszystkich kombinacji (isBest ≥ 0)? Symulacja pozwala rozstrzygnąć te pytania, bo dla każdego śladu znane jest jego pochodzenie.');

{
  const s = cs('Podsumowanie', 'W skrócie: wzór na tło działa dla wszystkich kombinacji, a przy isBest == 1 zawodzi',
    'Cztery liczby u góry dotyczą przedziału masy od 0,14 do 0,80 GeV/c². Przy wszystkich kombinacjach wzór daje 0,938 tła prawdziwego, a sygnał jest o 2,1 % za duży. Przy isBest == 1 wzór daje 1,207 tła prawdziwego, a sygnał jest o 5,5 % za mały. W przedziale 0,30–0,56 GeV/c² błąd sygnału sięga 12 %. Przyczyna: wybór jednej pary usuwa 23 % tła prawdziwego e⁺e⁻, a wartość wzoru zmienia się o 0,5 %, bo pary jednoimienne zwykle nie mają w zdarzeniu rywala. W zdarzeniach z rywalizacją wybór nie wyróżnia pary prawdziwej: zostaje jej tyle samo co tła. Użycie wszystkich kombinacji daje poprawny sygnał, a kosztuje tylko 2 % istotności statystycznej. Symulacja służy tu wyłącznie do zbadania wpływu isBest, bo tła kombinatorycznego w danych nie odtwarza. Część czwarta wyznacza poprawkę k z mieszania zdarzeń w danych, a część piąta pokazuje te same wielkości w danych eksperymentalnych.');
  const w = 2.9, g = (CW - 4*w)/3, y = 1.55;
  stat(s, '0,938 ± 0,004', 'wzór / tło prawdziwe, wszystkie kombinacje (isBest ≥ 0)', MX, y, w, 1.3, A);
  stat(s, '1,207 ± 0,006', 'wzór / tło prawdziwe, jedna para na zdarzenie (isBest == 1)', MX + (w + g), y, w, 1.3, B_);
  stat(s, '+2,1 %', 'błąd sygnału „wszystkie − wzór”, isBest ≥ 0', MX + 2*(w + g), y, w, 1.3, A);
  stat(s, '−5,5 %', 'błąd sygnału „wszystkie − wzór”, isBest == 1', MX + 3*(w + g), y, w, 1.3, B_);
  cap(s, `Masa pary M_{ee} od 0,14 do 0,80 ${M}, bez cięcia na kąt otwarcia, wszystkie zdarzenia.`, MX, 2.9, CW, 0.28);
  bl(s, [
    '**Przyczyna błędu.** Hipotezy e⁺e⁻ rywalizują w zdarzeniu i znika 23 % tła prawdziwego. Para jednoimienna zwykle nie ma rywala, więc wzór traci tylko 0,5 %.',
    '**isBest == 1 to para o najlepiej zgodnych czasach, a nie para najlepsza fizycznie.** Gdy ranking wybierze parę tła, para prawdziwa wypada z próbki i odjęcie tła jej nie przywróci: w danych w 0,14–0,80 GeV/c² tracimy 2–3 % sygnału.',
    `**Skutek w symulacji.** Zwykły wzór przy isBest == 1 odejmuje za dużo tła, więc sygnał jest zaniżony: o 5,5 % w 0,14–0,80 ${M} i o 12 % w 0,30–0,56 ${M}.`,
    `**Poprawka k z danych.** Także przy isBest ≥ 0 wzór jest za mały poniżej 0,30 ${M}. Z mieszania zdarzeń w danych k = ${T.k_int.k} ± ${T.k_int.err} w 0,14–0,80 ${M}.`,
    `**Dane, 0,14–0,80 ${M}.** Zwykły wzór przy isBest == 1 odejmuje o ${cbTooBig} % za dużo tła. Tła jest tu ponad dwa razy więcej niż sygnału, więc sygnał wychodzi o ${num2(D4.SbN/D4.SkN)} % za mały. To błąd wzoru, a nie utrata par: z poprawnym tłem różnica to tylko 2 %.`,
  ], { x: MX, y: 3.22, w: CW, h: 2.8, fontSize: 14.5, paraSpaceAfter: 5 });
  result(s, 'widmo inkluzywne: wszystkie kombinacje; isBest == 1 wymaga poprawki tła i traci część sygnału.');
}

// ================================================================ część 1
divider('Część 1: próbka', 'Część 1: symulacja, kanały i identyfikacja leptonów',
  ['Próbka i jej statystyka', 'Kanały rozpadu i pochodzenie π⁰', 'Identyfikacja leptonów w zmiennych β_{new} i q·p', 'Komplet cięć'],
  'W pierwszej części opisano, z czego zbudowana jest próbka: jakie są źródła par e⁺e⁻, ile ich jest i jak wybrano leptony.');

{
  const s = cs(P1('Próbka'), 'Próbka: 8,04 mln par e⁺e⁻ z symulacji SMASH i pełna historia każdego śladu',
    'Plik smash_std_sim_epem.root zawiera produkcję standardową SMASH, bez wag: liczba par jest bezpośrednio liczbą zrekonstruowanych par, żaden rozpad nie jest sztucznie wzmocniony. Po przejściu przez GEANT i rekonstrukcję zapisano trzy drzewa par: EpEm_ID, EpEp_ID i EmEm_ID. Każdy ślad ma w ntupli rodzaj cząstki, rodzica, proces i punkt powstania, dlatego można bezpośrednio policzyć, które pary są prawdziwe, a które należą do tła. Liczby na slajdzie są po cięciu na wierzchołek. Liczba zdarzeń wygenerowanych nie jest zapisana w tym pliku, dlatego statystyka podana jest wyłącznie dla zdarzeń zrekonstruowanych. Zdarzeń z parą jednoimienną, ale bez pary e⁺e⁻ jest 172 390: to głównie zdarzenia z dwoma leptonami tego samego znaku.');
  bl(s, [
    '**Generator:** SMASH, produkcja standardowa (gen4), bez wag, więc żaden rozpad nie jest sztucznie wzmocniony.',
    '**Detektor i rekonstrukcja:** GEANT, potem PAT; zmienna β_{new} i osobne cięcia PID dla symulacji.',
    '**Zapis:** trzy rodzaje par, e⁺e⁻, e⁺e⁺ oraz e⁻e⁻, bez powtórzeń kombinacji.',
    '**Historia śladu:** dla każdego śladu znane są rodzaj cząstki, rodzic, proces i punkt powstania.',
    'Liczby dotyczą par zrekonstruowanych. Liczby zdarzeń wygenerowanych nie ma w pliku.',
    '**Rola symulacji:** zbadanie wpływu isBest. Tła kombinatorycznego w danych symulacja nie odtwarza.',
  ], { x: MX, y: 1.55, w: 6.5, h: 4.4, fontSize: 15 });
  tbl(s, [
    ['Rodzaj pary', { t: 'Zdarzenia', ...num }, { t: 'Pary', ...num }],
    ['e⁺e⁻', { t: '7 951 064', ...num }, { t: '8 039 114', ...num }],
    ['e⁺e⁺', { t: '62 978', ...num }, { t: '63 132', ...num }],
    ['e⁻e⁻', { t: '192 886', ...num }, { t: '194 268', ...num }],
  ], { x: 7.55, y: 1.6, w: 5.18, colW: [1.78, 1.7, 1.7], rowH: 0.42, fontSize: 16 });
  cap(s, 'Po cięciu na wierzchołek, wszystkie zdarzenia, bez cięcia na kąt otwarcia. Zdarzenia z parą jednoimienną, ale bez pary e⁺e⁻: 172 390.', 7.55, 3.45, 5.18, 0.7);
  result(s, 'pary jednoimienne to tylko 3 % wszystkich par, a właśnie z nich wzór wylicza tło.');
}

{
  const s = cs(P1('Kanały'), 'Kanały: pary prawdziwe pochodzą z rozpadów Dalitza π⁰ i η oraz z konwersji fotonów',
    'Para prawdziwa to dwa leptony o przeciwnych znakach z tego samego punktu, z tym samym rodzicem i procesem. Konwersje fotonów też spełniają tę definicję. Rozpady Dalitza π⁰ dają 3,38 mln par, Dalitza η 2,01 mln, konwersje fotonów łącznie 2,40 mln. Konwersje mają bardzo małe kąty otwarcia, więc cięcie OA > 9° zostawia ich zaledwie 2718. Rozpad Dalitza π⁰ kończy się przy masie π⁰, powyżej ciągnie się Dalitz η do masy około 0,55 GeV/c². W tej próbce nie ma par e⁺e⁻ z ω ani z ρ, dlatego powyżej 0,56 GeV/c² nie ma par prawdziwych, zostaje samo tło.');
  tbl(s, [
    ['Kanał', 'Rozpad', { t: 'Bez\ncięcia', ...num }, { t: 'OA > 2°', ...num }, { t: 'OA > 4°', ...num }, { t: 'OA > 9°', ...num }],
    ['Dalitz π⁰', 'π⁰ → γ e⁺e⁻', { t: '3 376 070', ...num }, { t: '2 473 188', ...num }, { t: '1 432 703', ...num }, { t: '595 506', ...num }],
    ['Dalitz η', 'η → γ e⁺e⁻', { t: '2 005 578', ...num }, { t: '1 580 818', ...num }, { t: '1 134 115', ...num }, { t: '746 781', ...num }],
    ['konwersja γ z π⁰', 'γ → e⁺e⁻', { t: '2 287 253', ...num }, { t: '1 057 122', ...num }, { t: '137 928', ...num }, { t: '2 627', ...num }],
    ['konwersja γ z η', 'γ → e⁺e⁻', { t: '105 523', ...num }, { t: '44 036', ...num }, { t: '4 915', ...num }, { t: '84', ...num }],
    ['konwersja, inne γ', 'γ → e⁺e⁻', { t: '8 156', ...num }, { t: '3 366', ...num }, { t: '377', ...num }, { t: '7', ...num }],
    ['Razem', '', { t: '7 782 580', ...num }, { t: '5 158 530', ...num }, { t: '2 710 038', ...num }, { t: '1 345 005', ...num }],
  ], { x: MX, y: 1.55, w: 5.95, colW: [1.4, 0.95, 0.9, 0.9, 0.9, 0.9], rowH: 0.33, fontSize: 12, boldRows: [6], ruleRows: [5], margin: [2, 4, 2, 4] });
  cap(s, 'Pary e⁺e⁻, isBest ≥ 0, M_{ee} do 1,40 GeV/c^{2}, po cięciu na wierzchołek.', MX, 4.17, 5.95, 0.28);
  bl(s, [
    '**Rysunek:** ile par prawdziwych z każdego źródła przypada na jednostkę masy.',
    'Konwersje też są parami prawdziwymi: oba leptony pochodzą z tego samego fotonu.',
    `W tej próbce nie ma e⁺e⁻ z ω ani ρ, więc powyżej 0,56 ${M} zostaje samo tło.`,
  ], { x: MX, y: 4.44, w: 5.95, h: 1.58, fontSize: 13.5, paraSpaceAfter: 3 });
  fig(s, 'figs/channels.png', 6.75, 1.5, 5.98, 4.5, 'Rozkłady masy par prawdziwych w czterech grupach kanałów, skala logarytmiczna, punkty połączone cienkimi liniami. Dalitz η sięga najwyższej masy, Dalitz π⁰ kończy się przy masie π⁰, konwersje leżą przy najmniejszych masach.');
  result(s, 'konwersje to 31 % par bez cięcia, ale po OA > 9° praktycznie znikają.');
}

{
  const s = cs(P1('Pochodzenie π⁰'), 'Pochodzenie π⁰ w rozpadach Dalitza: około 90 % zapisano jako pierwotne, około 6 % pochodzi z ρ i ω',
    'Kod geninfo mówi, skąd pochodzi przodek leptonu: dwukrotnie powtórzony numer cząstki oznacza cząstkę pierwotną, a długi kod oznacza rozpad zapisany jako łańcuch, na przykład 420807 to ρ⁺ → π⁺π⁰. Kod 52080908 jest zaokrąglonym 52080907: kody są przechowywane w zmiennej typu float, która dla liczb powyżej 2 do potęgi 24 traci dokładność. 89,9 % π⁰ z rozpadów Dalitza zapisano jako pierwotne. Rozpady ρ i ω dają razem około 6 %. Dla η odsetek zapisanych jako pierwotne to 99,99 %. Dla konwersji fotonów rozkład źródeł jest bardzo podobny, z nieco większym udziałem η.');
  const rows = [['Kod geninfo przodka', { t: `Dalitz π⁰ (${T.src_tot.pi0}) [%]`, ...num }, { t: `Konwersje γ (${T.src_tot.conv}) [%]`, ...num }]];
  T.src.forEach(r => rows.push([r.lab, { t: r.pi0, ...num }, { t: r.conv, ...num }]));
  tbl(s, rows, { x: MX, y: 1.6, w: 8.1, colW: [3.5, 2.3, 2.3], rowH: 0.45, fontSize: 15, ruleRows: [5] });
  cap(s, 'Pary prawdziwe e⁺e⁻, isBest ≥ 0, bez cięcia na kąt otwarcia. Kolumna Dalitz π⁰: skąd pochodzi π⁰. Kolumna konwersji: mezon, z którego pochodzi foton.', MX, 5.05, 8.1, 0.5);
  stat(s, '99,99 %', 'par z Dalitza η ma η zapisane jako pierwotne', 9.0, 1.6, 3.73, 1.35);
  stat(s, 'około 6 %', 'π⁰ z rozpadów Dalitza pochodzi z rozpadów ρ i ω', 9.0, 3.1, 3.73, 1.35);
  result(s, 'około 90 % π⁰ z rozpadów Dalitza i 99,99 % η zapisano jako cząstki pierwotne.');
}

{
  const s = cs(P1('Identyfikacja leptonów'), 'Identyfikacja leptonów: β_{new} względem q·p, cięcie 5σ',
    'Oś pozioma to ładunek razy pęd, więc elektrony leżą po lewej, a pozytony po prawej stronie; kontury dla e⁻ są odbite symetrycznie na lewą stronę. Kontury 1σ, 3σ i 5σ pochodzą z pliku pp45_pid_cuts_sim_betanew_v6.root, wyznaczonego osobno dla symulacji tym samym algorytmem co w danych. W wersji v6 kontur e⁺ jest praktycznie konturem e⁻, bo kontur e⁺ fluktuował. Zastosowane cięcie to 5σ, nałożone na każdy z dwóch leptonów pary już przy zapisie pliku. Liczby na slajdach nie były przeliczane po zmianie konturów, bo różnice są promilowe: histogram i procenty w tabeli pochodzą z produkcji z konturami v5, dlatego pojedyncze ślady mogą wystawać poza kontur v6. W 3σ mieści się 98,8 % śladów, w 1σ około 70 %. Histogram zawiera każdy ślad raz: 7 980 135 e⁺ i 8 005 683 e⁻.');
  fig(s, 'figs/pid.png', 6.5, 1.55, 6.23, 4.42, 'Dwuwymiarowy rozkład beta_new względem iloczynu ładunku i pędu. Elektrony po lewej, pozytony po prawej. Paleta od niebieskiej przy małej liczbie śladów do żółtej i czerwonej w maksimum. Czarne kontury z pliku v6: 1 sigma kropkami, 3 sigma kreskami, 5 sigma linią ciągłą.');
  bl(s, [
    '**Rysunek:** każdy punkt to ślad leptonu, czarne linie to granice cięć.',
    '**β_{new}:** prędkość wyznaczona z czasu względnego.',
    '**q·p:** ładunek razy pęd, więc e⁻ leżą po lewej, a e⁺ po prawej stronie.',
    '**Kontury v6:** osobne dla symulacji, algorytm jak w danych; kontur e⁺ praktycznie równy e⁻, odbity symetrycznie.',
    '**Cięcie graficzne:** 5σ na każdy z dwóch leptonów, już przy zapisie.',
  ], { x: MX, y: 1.55, w: 5.6, h: 2.75, fontSize: 14, paraSpaceAfter: 4 });
  tbl(s, [
    ['Kontur', { t: 'e⁺', ...num }, { t: 'e⁻', ...num }],
    ['1σ', { t: '70,4 %', ...num }, { t: '69,9 %', ...num }],
    ['3σ', { t: '98,8 %', ...num }, { t: '98,8 %', ...num }],
    ['5σ', { t: '100 %', ...num }, { t: '100 %', ...num }],
  ], { x: MX, y: 4.3, w: 3.6, colW: [1.2, 1.2, 1.2], rowH: 0.32, fontSize: 14, boldRows: [3] });
  cap(s, 'Ślady z par e⁺e⁻, każdy raz. Procenty z konturami v5, różnice względem v6 są promilowe.', 4.4, 4.35, 1.8, 1.5);
  result(s, 'w 1σ mieści się około 70 % śladów, w 3σ 98,8 %, analiza używa cięcia 5σ.');
}

{
  const s = cs(P1('Selekcja'), 'Selekcja: wszystkie użyte cięcia',
    'Cięcie na wierzchołek jest pierwszym krokiem selekcji i jest konieczne: reakcje mają pochodzić z tarczy, a brak zrekonstruowanego wierzchołka (wartość −1000) oznacza odrzucenie. W tej symulacji usuwa ono 3,9 % par e⁺e⁻, 3,0 % e⁺e⁺ i 5,2 % e⁻e⁻. Cięcie PID 5σ jest już w pliku. Porównywane są dwa warianty isBest. Główne liczby są bez cięcia na kąt otwarcia; cięcia OA > 2°, 4° i 9° sprawdzono osobno. Wzór na tło liczony jest w każdym przedziale masy osobno, nigdy z sum po szerokich przedziałach.');
  tbl(s, [
    ['Cięcie', 'Wartość', 'Gdzie', 'Po co', 'Co zmienia'],
    ['Wierzchołek', 'eVertReco_z > −500 mm', 'zdarzenie, jako pierwsze', 'reakcje z tarczy; −1000 oznacza brak wierzchołka', 'usuwa 3,9 % par e⁺e⁻, 3,0 % e⁺e⁺, 5,2 % e⁻e⁻'],
    ['PID', 'wewnątrz konturu 5σ w (q·p, β_{new})', 'każdy lepton pary, przy zapisie', 'wybór leptonów', 'w pliku są tylko ślady wewnątrz 5σ'],
    ['Wybór pary', 'isBest ≥ 0 albo isBest == 1', 'każda para, osobno dla e⁺e⁻, e⁺e⁺, e⁻e⁻', 'wszystkie kombinacje albo jedna para na zdarzenie', 'to główny temat części 3'],
    ['Kąt otwarcia', 'brak; osobno OA > 2°, 4°, 9°', 'każda para', 'usuwa konwersje i pary ze wspólnym pierścieniem RICH', 'w liczbach głównych nie ma cięcia OA'],
  ], { x: MX, y: 1.55, w: CW, colW: [1.5, 2.95, 2.35, 2.95, 2.38], rowH: [0.38, 0.78, 0.78, 0.78, 0.78], fontSize: 14 });
  cap(s, `Masa pary: 58 przedziałów od 0 do 1,40 ${M} (co 2,5 MeV/c^{2} do 0,05 ${M}, co 20 MeV/c^{2} od 0,14 do 0,80 ${M}). Wszystkie zdarzenia, bez wyboru trygera. Wszystkie liczby na slajdach są po cięciu na wierzchołek i po cięciu PID 5σ.`, MX, 5.3, CW, 0.7);
}

// ================================================================ część 2
divider('Część 2: sygnał i tło', 'Część 2: sygnał i tło przy wszystkich kombinacjach',
  ['Definicje, tło prawdziwe i warunki poprawności wzoru', 'Widmo e⁺e⁻ i stosunek sygnału do tła, także w skali liniowej', 'Porównanie wzoru z tłem policzonym ze śladów i przyczyny różnicy poniżej 0,30 GeV/c^{2}'],
  'W drugiej części używane są wszystkie zapisane kombinacje, czyli isBest ≥ 0. To punkt odniesienia dla części trzeciej.');

{
  const s = cs(P2('Definicje'), 'Cztery wielkości na widmie masy e⁺e⁻ i sposób ich liczenia',
    'N+− to wszystkie pary e⁺e⁻. Sygnał prawdziwy i tło prawdziwe liczone są z historii śladów: sygnał to pary z jednego punktu, jednego rodzica i jednego procesu, a tło to pary, których ślady mają różne źródła, oraz pary z leptonem wtórnym albo z fałszywym śladem. Tło ze wzoru liczone jest z par jednoimiennych, w każdym przedziale masy osobno, bez poprawki akceptancji, czyli z k = 1. Osobną grupą jest jedna cząstka zrekonstruowana dwa razy: 40 907 par poniżej 0,14 GeV/c², które nie należą ani do sygnału, ani do tła prawdziwego i dlatego zostają w N − CB.');
  const w = 5.95, g = CW - 2*w, h = 1.72, y1 = 1.55, y2 = y1 + h + 0.15;
  card(s, 'N_{+−}: wszystkie pary e⁺e⁻', 'Wszystkie zrekonstruowane pary w przedziale masy, bez użycia historii śladów.', MX, y1, w, h);
  card(s, 'Sygnał prawdziwy', 'Oba ślady to leptony o przeciwnych znakach, z tego samego punktu, z tym samym rodzicem i procesem. Liczony z historii śladów.', MX + w + g, y1, w, h);
  card(s, 'Tło prawdziwe', 'Pozostałe pary: leptony z różnych źródeł, lepton wtórny albo fałszywy ślad. Liczone z historii śladów.', MX, y2, w, h);
  card(s, 'Tło ze wzoru', 'Z par jednoimiennych, w każdym przedziale masy osobno, bez poprawki akceptancji (k = 1).', MX + w + g, y2, w, h);
  eq(s, 'cb', MX + w + g + 0.15, y2 + 1.17, 1.9, 'Tło ze wzoru równa się dwa razy pierwiastek z iloczynu liczb par e⁺e⁺ i e⁻e⁻.');
  eq(s, 's', MX + w + g + 3.1, y2 + 1.2, 1.9, 'Sygnał ze wzoru równa się liczba wszystkich par e⁺e⁻ minus tło ze wzoru.');
  tx(s, `Osobna grupa: jedna cząstka zrekonstruowana dwa razy, 40 907 par poniżej 0,14 ${M}. Nie jest ani sygnałem, ani tłem prawdziwym, więc zostaje w różnicy N_{+−} − CB.`, { x: MX, y: 5.3, w: CW, h: 0.65 });
}

{
  const s = cs(P2('Tło prawdziwe i warunki wzoru'), 'Tło prawdziwe ze śladów i warunki poprawności wzoru na tło',
    'Tło prawdziwe liczone jest para po parze z historii śladów. Dla każdego śladu symulacja podaje rodzaj cząstki, rodzica, proces i punkt powstania. Para jest sygnałem, gdy oba ślady to lepton i antylepton z tego samego punktu, z tym samym rodzicem i procesem. Tłem jest każda inna para: leptony z różnych źródeł, para z leptonem wtórnym i para z fałszywym śladem. Jedna cząstka zrekonstruowana dwa razy to osobna grupa. Wzór wynika z prostego rachunku. Jeśli w zdarzeniu jest N niezależnych źródeł par i każde daje jeden e⁺ i jeden e⁻, a prawdopodobieństwa rejestracji to ε₊ i ε₋, to liczba par e⁺e⁻ z różnych źródeł wynosi ⟨N(N−1)⟩ε₊ε₋, a liczby par jednoimiennych ⟨N(N−1)⟩ε₊²/2 oraz ⟨N(N−1)⟩ε₋²/2. Stąd tożsamość, niezależna od rozkładu N. Działa ona, gdy spełnione są trzy warunki z tabeli. Warunek pierwszy jest złamany przez elektrony wtórne, które powstają tylko jako e⁻, i przez fałszywe ślady, częstsze dla ładunku dodatniego. Warunek drugi wymaga, żeby akceptancja par różnoimiennych była średnią geometryczną akceptancji par jednoimiennych; sam stały stosunek N₊₊/N₋₋ różny od 1 mu nie przeszkadza, przeszkadza inny kształt akceptancji e⁺ i e⁻ w funkcji masy. Jego złamanie opisuje poprawka k. Warunek trzeci łamie wybór isBest == 1.');
  const lw = 5.75, rx = MX + lw + 0.25, rw = CW - lw - 0.25;
  card(s, 'Tło prawdziwe: para po parze, z historii śladów', [
    '**Dane o śladzie:** rodzaj cząstki, rodzic, proces, punkt powstania.',
    '**Sygnał:** e⁺ i e⁻ z tego samego punktu, z tym samym rodzicem i procesem.',
    '**Tło:** leptony z różnych źródeł (Dalitz + Dalitz, Dalitz + konwersja, konwersja + konwersja), para z leptonem wtórnym albo z fałszywym śladem.',
    '**Osobno:** jedna cząstka zrekonstruowana dwa razy.',
    '**Liczenie:** w każdym przedziale masy, osobno dla isBest ≥ 0 i isBest == 1.',
  ], MX, 1.55, lw, 4.4, { size: 15, gap: 8, fill: C.accent1 });
  tx(s, '**Skąd wzór:** N niezależnych źródeł w zdarzeniu, każde daje jeden e⁺ i jeden e⁻; ε_{+}, ε_{−}: prawdopodobieństwa rejestracji', { x: rx, y: 1.55, w: rw, h: 0.55, fontSize: 14 });
  eq(s, 'derive', rx, 2.14, 1.42, 'N plus plus równa się połowa średniej N razy N minus 1 razy epsilon plus do kwadratu, analogicznie N minus minus, a tło par różnoimiennych równa się średnia N razy N minus 1 razy epsilon plus razy epsilon minus, czyli dwa pierwiastki z iloczynu N plus plus i N minus minus.');
  tbl(s, [
    ['Warunek', 'W tej próbce'],
    ['1. Każde źródło daje tyle samo e⁺ co e⁻', 'Dalitz i konwersje: tak. Elektrony wtórne (tylko e⁻) i fałszywe ślady (głównie e⁺): nie'],
    ['2. Akceptancja: A_{+−} = √(A_{++} A_{−−}) w każdym przedziale masy', `Nie poniżej 0,30 ${M}: N_{++}/N_{−−} rośnie z ${T.why[1].ppmm} do ${T.why[3].ppmm}; opisuje to poprawka k`],
    ['3. Ta sama selekcja par +− i par jednoimiennych', 'isBest ≥ 0: tak. isBest == 1: nie (część 3)'],
  ], { x: rx, y: 3.32, w: rw, colW: [2.95, rw - 2.95], rowH: [0.34, 0.72, 0.78, 0.5], fontSize: 13 });
  result(s, 'wzór jest ścisły przy symetrii ładunku źródeł, zgodnej akceptancji i tej samej selekcji par.');
}

{
  const s = cs(P2('Widmo'), 'Widmo e⁺e⁻ przy isBest ≥ 0: tło ze wzoru odejmuje się niemal tak dobrze jak tło policzone ze śladów',
    'Wszystkie krzywe na tym slajdzie i na dwóch następnych są dla isBest ≥ 0, czyli dla wszystkich zapisanych kombinacji. Czarne punkty to wszystkie pary e⁺e⁻. Niebieska linia to sygnał prawdziwy, otwarte kółka to różnica wszystkich par i tła prawdziwego, fioletowe trójkąty to różnica wszystkich par i tła ze wzoru. Te trzy zbiory leżą na sobie aż do około 0,5 GeV/c². W dolnym pasku pokazano stosunek sygnału do tła: wypełniony histogram dla tła prawdziwego, trójkąty dla wzoru. Sygnał przeważa nad tłem do około 0,43 GeV/c². Powyżej 0,52 GeV/c² sygnału praktycznie nie ma, a w 0,56–0,80 GeV/c² jest już tylko tło: 1780 par, z czego prawdziwych 5. Wzór daje tam 1783 pary tła wobec 1775 prawdziwych.');
  fig(s, 'figs/spectrum.png', MX, 1.5, 6.5, 4.5, 'Widmo masy e⁺e⁻, skala logarytmiczna: wszystkie pary, sygnał prawdziwy, wszystkie pary minus tło prawdziwe, wszystkie pary minus tło ze wzoru. Pod spodem pasek ze stosunkiem sygnału do tła.', { left: true });
  bl(s, [
    '**Selekcja:** isBest ≥ 0, czyli wszystkie kombinacje.',
    '**Czarne punkty:** wszystkie pary e⁺e⁻, N_{+−}.',
    '**Niebieska linia:** sygnał prawdziwy, z historii śladów.',
    '**Kółka:** N_{+−} minus tło prawdziwe.',
    '**Trójkąty:** N_{+−} minus tło ze wzoru.',
    `**Dolny pasek:** sygnał / tło. Sygnał przeważa do około 0,43 ${M}.`,
  ], { x: 7.3, y: 1.55, w: 5.43, h: 3.65, fontSize: 15 });
  cap(s, 'Bez cięcia na kąt otwarcia, wszystkie zdarzenia. Przedziały masy różnej szerokości, dlatego liczby par podzielono przez szerokość przedziału.', 7.3, 5.2, 5.43, 0.72);
  result(s, `sygnał ze wzoru zgadza się z prawdziwym w 2,1 % w przedziale 0,14–0,80 ${M}.`);
}

{
  const s = cs(P2('Skala liniowa'), `To samo widmo w skali liniowej powyżej 0,1 ${M}: różnicę tła widać dopiero tutaj`,
    'To dokładnie te same dane co na poprzednim slajdzie, isBest ≥ 0, w skali liniowej i tylko powyżej 0,1 GeV/c². Dodano tło prawdziwe (czerwona linia) i tło ze wzoru (czerwone kwadraty), żeby zobaczyć ich różnicę bezpośrednio. W przedziale 0,14–0,30 GeV/c² kwadraty leżą pod linią: wzór jest za mały o około 9 %. Różnica ta jest mała względem sygnału, dlatego w skali logarytmicznej jej nie widać. Sygnał ze wzoru leży przez to minimalnie powyżej sygnału prawdziwego. W dolnym pasku pokazano stosunek sygnału do tła w skali liniowej.');
  fig(s, 'figs/spectrum_lin.png', MX, 1.5, 8.0, 4.5, 'Widmo masy e⁺e⁻ w skali liniowej powyżej 0,1 GeV/c²: wszystkie pary, sygnał prawdziwy, tło prawdziwe, tło ze wzoru i dwa warianty sygnału. Pod spodem stosunek sygnału do tła w skali liniowej.', { left: true });
  bl(s, [
    '**Selekcja:** isBest ≥ 0, jak na poprzednim slajdzie.',
    `**Czerwone:** tło prawdziwe (linia) i ze wzoru (kwadraty). W 0,14–0,30 ${M} kwadraty leżą pod linią.`,
    '**Niebieskie i fioletowe:** sygnał prawdziwy i ze wzoru; różnica jest mała, bo sygnał jest tu większy od tła.',
    `**Pasek:** sygnał / tło spada poniżej 1 około 0,43 ${M}.`,
  ], { x: 8.8, y: 1.55, w: 3.93, h: 4.4, fontSize: 15 });
  result(s, `skala liniowa pokazuje, że tło ze wzoru jest za małe w 0,14–0,30 ${M}.`);
}

{
  const s = cs(P2('Wzór a tło ze śladów'), `To samo widmo w liczbach, isBest ≥ 0: powyżej 0,30 wzór się zgadza, w 0,14–0,30 ${M} jest za mały o 9 %`,
    'To ta sama próbka co na dwóch poprzednich slajdach, czyli isBest ≥ 0; zmienia się tylko sposób pokazania: tam widma, tu ich stosunek. Skala logarytmiczna nie pokazuje różnicy 9 %, bo sygnał jest w tym przedziale prawie 4 razy większy od tła. Stosunek tła ze wzoru do tła prawdziwego jest bliski jedności w szerokim zakresie mas. Poniżej 0,02 GeV/c² wzór jest za duży o 23 %, a w tym obszarze występują pary z fałszywym śladem: poniżej 0,14 GeV/c² 70 % z nich dzieli pierścień RICH z partnerem. Przy cięciu na kąt otwarcia stosunek spada: 1,11 dla OA > 2° i 1,03 ± 0,04 dla OA > 4°. Mimo tego błędu sygnał w tym przedziale zmienia się o mniej niż 1 %, bo stosunek sygnału do tła wynosi tam około 390. W przedziale 0,14–0,30 GeV/c² wzór jest o 9 % za mały. Z definicji jest to poprawka akceptancji k, czyli stosunek tła prawdziwego do wzoru. W części drugiej i trzeciej przyjęto k = 1, żeby pokazać sam wzór; poprawkę wyznaczono w części czwartej.');
  const rows = [['Przedział [GeV/c^{2}]', { t: 'Tło\nprawdziwe', ...num }, { t: 'Wzór', ...num }, { t: 'Wzór /\nprawdziwe', ...num }]];
  T.windows.forEach(r => rows.push([r.w, { t: r.all.B, ...num }, { t: r.all.CB, ...num }, { t: `${r.all.R} ± ${r.all.eR}`, ...num }]));
  tbl(s, rows, { x: MX, y: 1.55, w: 6.2, colW: [1.9, 1.3, 1.1, 1.9], rowH: 0.36, fontSize: 14, boldRows: [6], ruleRows: [5] });
  bl(s, [
    '**Ta sama próbka co wyżej:** isBest ≥ 0. W skali logarytmicznej 9 % nie widać, bo sygnał jest prawie 4 razy większy od tła.',
    `Poniżej 0,02 ${M} wzór jest za duży o 23 %; przy OA > 4° stosunek wynosi 1,03 ± 0,04.`,
    `W 0,14–0,30 ${M} wzór jest za mały o 9 %, przyczyna na następnym slajdzie.`,
  ], { x: MX, y: 4.5, w: 6.2, h: 1.5, fontSize: 14 });
  fig(s, 'figs/ratio_all.png', 7.05, 1.55, 5.68, 4.4, 'Stosunek tła ze wzoru do tła prawdziwego w funkcji masy pary dla wszystkich kombinacji. Wartości blisko jedności, minimum około 0,85 przy 0,15 GeV/c², na początku widma powyżej 1.');
  cap(s, 'Rysunek: ile razy tło ze wzoru różni się od tła prawdziwego; 1 oznacza zgodność.', 7.05, 5.55, 5.68, 0.4);
  result(s, `w przedziale 0,14–0,80 ${M} wzór daje 0,938 ± 0,004 tła prawdziwego.`);
}

{
  const s = cs(P2('Przyczyna różnicy'), `Dlaczego poniżej 0,30 ${M} wzór jest za mały: warunek akceptancji nie jest spełniony`,
    'Tabela sprawdza warunki ze slajdu 11 osobno dla grup źródeł. Stosunek wzoru do tła dla samych leptonów z różnych źródeł, czyli dla par, które formalnie spełniają warunek pierwszy, jest prawie taki sam jak dla całego tła: 0,874 wobec 0,872 w 0,14–0,20 GeV/c². Deficyt nie pochodzi więc od elektronów wtórnych ani od fałszywych śladów; ich udział w N₋₋ częściowo równoważy ich brak w modelu. Mieszanie zdarzeń, w którym nie ma par skorelowanych, również daje k większe od 1. To pokazuje, że przyczyną jest akceptancja, czyli warunek drugi. Mieszanie wyjaśnia od jednej trzeciej do połowy różnicy, zależnie od puli leptonów; pula ze zdarzeń z co najmniej trzema leptonami daje większą poprawkę niż pula ze wszystkich zdarzeń. Wynika to z tego, że pula ze wszystkich zdarzeń składa się głównie z leptonów par prawdziwych. Pełną poprawkę wyznaczono w części czwartej. Powyżej 0,30 GeV/c² wzór zgadza się w granicach 1 %.');
  const rows = [['Przedział\n[GeV/c²]', { t: 'Wzór /\nprawdziwe', ...num }, { t: 'Różne źródła:\nwzór / prawdz.', ...num }, { t: 'e⁺e⁺ /\ne⁻e⁻', ...num }, { t: 'k_{mix}\nwszystkie', ...num }, { t: 'k_{mix}\nn ≥ 3', ...num }]];
  T.why.forEach(r => rows.push([r.w, { t: r.R, ...num, bold: true }, { t: r.Rlep, ...num }, { t: r.ppmm, ...num }, { t: r.kA, ...num }, { t: r.kN, ...num }]));
  tbl(s, rows, { x: MX, y: 1.55, w: 6.75, colW: [1.15, 1.15, 1.45, 0.85, 1.1, 1.05], rowH: [0.62, 0.36, 0.36, 0.36, 0.36], fontSize: 12.5 });
  bl(s, [
    '**To nie elektrony wtórne** (6–8 % tła): dla samych leptonów z różnych źródeł (Dalitz + Dalitz itd.) stosunek jest taki sam.',
    '**To akceptancja:** w mieszaniu nie ma par skorelowanych, a k_{mix} też jest większe od 1 (warunek 2).',
    '**Mieszanie to za mało:** tłumaczy od 1/3 do 1/2 różnicy, reszta zależy od tego, z jakich zdarzeń są leptony.',
  ], { x: MX, y: 3.92, w: 6.75, h: 2.1, fontSize: 15 });
  fig(s, 'figs/why_low.png', 7.55, 1.55, 5.18, 4.4, 'Stosunek tła ze wzoru do tła prawdziwego i odwrotność k z mieszania dla dwóch pul leptonów w funkcji masy. Poniżej 0,3 GeV/c² wszystkie krzywe leżą pod jednością, stosunek wzoru do tła najniżej.');
  cap(s, 'Rysunek: czarne punkty mówią, jak bardzo wzór różni się od prawdy; linie pokazują, ile z tego wyjaśnia mieszanie zdarzeń.', 7.55, 5.3, 5.18, 0.6);
  result(s, 'poniżej 0,30 wzór potrzebuje poprawki k; poniżej 0,14 i powyżej 0,30 wystarcza k bliskie 1.');
}

// ================================================================ część 3
divider('Część 3: isBest', 'Część 3: pary jednoimienne i wybór isBest',
  ['Skład par e⁺e⁺ i e⁻e⁻', 'Mechanizm działania isBest i klasy krotności leptonów', 'Dlaczego isBest == 1 psuje wzór na tło i sygnał', 'Do czego prowadzi isBest ≥ 0'],
  'W trzeciej części pokazano, skąd bierze się błąd wzoru przy wyborze jednej pary i do czego prowadzi użycie wszystkich kombinacji.');

{
  const s = cs(P3('Pary jednoimienne'), 'Pary e⁺e⁺ i e⁻e⁻ mają te same źródła co tło e⁺e⁻, lecz e⁻e⁻ ma dodatkowo elektrony wtórne',
    'Pary jednoimienne powstają głównie z połączenia leptonów z dwóch różnych rozpadów Dalitza lub konwersji. Proporcje są bardzo podobne jak w tle prawdziwym e⁺e⁻. Różnica dotyczy leptonów wtórnych: elektrony z rozproszenia Comptona i elektrony delta powstają wyłącznie jako e⁻. Dlatego w parach e⁻e⁻ udział leptonów wtórnych wynosi 13,1 %, a w e⁺e⁺ tylko 0,6 %. Tło e⁺e⁻ ma wartość pośrednią, 7,1 %, czego należy oczekiwać, bo wzór jest średnią geometryczną liczb par e⁺e⁺ i e⁻e⁻. Całkowita liczba par e⁻e⁻ jest trzy razy większa niż e⁺e⁺. Fałszywe ślady dominują w parach e⁺e⁺ przy najmniejszych masach: to głównie para ze wspólnym pierścieniem RICH.');
  fig(s, 'figs/like_sign.png', MX, 1.5, CW, 4.75, 'Rozkłady masy par e⁺e⁺ (po lewej) i e⁻e⁻ (po prawej) według źródeł: Dalitz plus Dalitz, Dalitz plus konwersja, konwersja plus konwersja, leptony wtórne i inne, fałszywy ślad.', { left: true });
  result(s, 'skład tła e⁺e⁻ leży między składem e⁺e⁺ a e⁻e⁻ (tabela w dodatku), więc wzór ma dobre podstawy.', 6.3, 0.45);
  cap(s, 'Rysunek: skąd pochodzą leptony w parach jednoimiennych; każdy kolor to inne źródło, oś pionowa jest logarytmiczna.', MX, 6.8, CW, 0.28);
}

function chip(s, label, sub, x, y, w, kept) {
  const head = rt(label, { bold: true, fontSize: 17, color: kept ? C.background1 : C.text2 });
  head[head.length - 1].options = { ...head[head.length - 1].options, breakLine: true };
  s.addText([...head, ...rt(sub, { fontSize: 12, color: kept ? C.background1 : C.text2 })],
    { shape: pres.ShapeType.roundRect, rectRadius: 0.06, x, y, w, h: 0.72, fill: { color: kept ? C.accent1 : C.background1 }, line: { color: kept ? C.accent1 : C.text2, width: 1, dashType: kept ? 'solid' : 'dash' }, align: 'center', valign: 'middle', margin: 0, objectName: kept ? 'Hipoteza zachowana' : 'Hipoteza odrzucona' });
}
{
  const s = cs(P3('Mechanizm'), 'Co robi isBest: w zdarzeniu rywalizują tylko hipotezy tego samego rodzaju',
    'Znacznik isBest wybiera w zdarzeniu jedną hipotezę: tę z najmniejszym χ² zgodności czasów przelotu. Rywalizacja odbywa się osobno w każdym rodzaju hipotez, a nie między rodzajami. Przykład po lewej: zdarzenie z dwoma pozytonami i jednym elektronem. Są dwie hipotezy e⁺e⁻, więc jedna odpada, ale jedyna hipoteza e⁺e⁺ nie ma rywala i zostaje. Przykład po prawej: dwa pozytony i dwa elektrony. Są cztery hipotezy e⁺e⁻, więc odpadają trzy, ale e⁺e⁺ i e⁻e⁻ są po jednej i zostają obie. Na tym polega istota problemu: pary jednoimienne nie mają konkurencji, która eliminowałaby je tak, jak eliminowane są pary e⁺e⁻.');
  const pw = 5.95, gap = CW - 2*pw, py = 1.55, ph = 3.3, xB = MX + pw + gap;
  card(s, 'Zdarzenie z 2 e⁺ i 1 e⁻', '', MX, py, pw, ph);
  card(s, 'Zdarzenie z 2 e⁺ i 2 e⁻', '', xB, py, pw, ph);
  const lab = (t, x, y) => tx(s, t, { x, y, w: pw - 0.4, h: 0.3, fontSize: 14, color: C.text2, name: 'Etykieta' });
  lab('Hipotezy e⁺e⁻: dwie, zostaje jedna', MX + 0.2, py + 0.55);
  chip(s, 'e⁺_{1} e⁻', 'isBest == 1', MX + 0.2, py + 0.88, 1.6, true);
  chip(s, 'e⁺_{2} e⁻', 'odrzucona', MX + 1.95, py + 0.88, 1.6, false);
  lab('Hipotezy e⁺e⁺: jedna, bez rywala', MX + 0.2, py + 1.95);
  chip(s, 'e⁺_{1} e⁺_{2}', 'isBest == 1', MX + 0.2, py + 2.28, 1.6, true);
  lab('Hipotezy e⁺e⁻: cztery, zostaje jedna', xB + 0.2, py + 0.55);
  [['e⁺_{1} e⁻_{1}', true], ['e⁺_{1} e⁻_{2}', false], ['e⁺_{2} e⁻_{1}', false], ['e⁺_{2} e⁻_{2}', false]].forEach(([t, k], i) => chip(s, t, k ? 'isBest == 1' : 'odrzucona', xB + 0.2 + i*1.4, py + 0.88, 1.27, k));
  lab('Hipotezy e⁺e⁺ i e⁻e⁻: po jednej, bez rywala', xB + 0.2, py + 1.95);
  chip(s, 'e⁺_{1} e⁺_{2}', 'isBest == 1', xB + 0.2, py + 2.28, 1.6, true);
  chip(s, 'e⁻_{1} e⁻_{2}', 'isBest == 1', xB + 1.95, py + 2.28, 1.6, true);
  tx(s, 'isBest == 1 ma hipoteza o najmniejszym χ² czasów przelotu wśród hipotez tego samego rodzaju w zdarzeniu. Przy isBest ≥ 0 zostają wszystkie hipotezy.', { x: MX, y: 5.05, w: CW, h: 0.8 });
  result(s, 'z dwóch lub czterech hipotez e⁺e⁻ zostaje jedna, a pary jednoimienne zostają wszystkie.');
}

{
  const s = cs(P3('Co zostaje'), 'isBest == 1 zachowuje prawie wszystkie pary prawdziwe i jednoimienne, ale tylko 77 % tła e⁺e⁻',
    'Wykres pokazuje, jaki ułamek par każdego rodzaju zostaje po isBest == 1. Zostaje 99,6 % par prawdziwych, ale nie dlatego, że isBest trafnie je wskazuje: 99 % z nich leży w zdarzeniach z dwoma leptonami, gdzie nie ma wyboru; następny slajd pokazuje zdarzenia z rywalizacją. Pary e⁺e⁺ i e⁻e⁻ zachowują się podobnie, ponad 99 % zostaje, bo w typowym zdarzeniu jest ich tylko jedna. Tło prawdziwe e⁺e⁻ traci znacznie więcej: w przedziale 0,14–0,80 GeV/c² zostaje 77,3 %, a zależnie od masy od około 73 % do 90 %. Wszystkie punkty dotyczą tych samych zdarzeń i tej samej selekcji.');
  fig(s, 'figs/kept_pm.png', MX, 1.5, 6.3, 4.55, 'Widma par e⁺e⁻, pary prawdziwe i tło prawdziwe, dla isBest ≥ 0 (linie) i isBest == 1 (punkty). Pod spodem ułamek par, które zostają: pary prawdziwe blisko 1, tło prawdziwe od 0,73 do 0,90.', { left: true });
  fig(s, 'figs/kept_ls.png', 7.15, 1.5, 5.58, 3.2, 'Widma par e⁺e⁺ i e⁻e⁻ dla isBest ≥ 0 (linie) i isBest == 1 (punkty). Pod spodem ułamek par, które zostają, blisko 1.', { left: true });
  bl(s, [
    '**Rysunki:** u góry ile jest par, na dole jaka ich część zostaje.',
    '**99,6 %** par prawdziwych zostaje: 99 % z nich nie ma rywala.',
    '**77,3 %** tła e⁺e⁻ zostaje, zależnie od masy 73–90 %.',
    '**Ponad 99 %** par e⁺e⁺ i e⁻e⁻ zostaje.',
  ], { x: 7.15, y: 4.78, w: 5.58, h: 1.25, fontSize: 13, paraSpaceAfter: 3 });
  result(s, 'isBest == 1 usuwa 23 % tła prawdziwego e⁺e⁻, a wartość wzoru zmienia o 0,5 %.');
}

{
  const s = cs(P3('Klasy krotności'), 'Strata tła e⁺e⁻ rośnie z liczbą leptonów: zostaje około 1 / (liczba hipotez e⁺e⁻)',
    'Klasa n to liczba zidentyfikowanych leptonów w zdarzeniu, czyli suma e⁺ i e⁻. Przy n = 2 jest jedna hipoteza e⁺e⁻ i nie ma rywalizacji. Przy n = 3, układ 2 + 1, są dwie hipotezy e⁺e⁻: zostaje jedna, czyli około 0,5. Przy n = 4, układ 2 + 2, są cztery hipotezy e⁺e⁻: zostaje około 0,25. Pary jednoimienne zostają prawie wszystkie, bo w tych układach są pojedyncze. Klasa n = 2 daje 56 % tła prawdziwego i nie traci nic. To znaczy, że wzrost stosunku wzoru do tła pochodzi wyłącznie z pozostałych 44 % tła, czyli ze zdarzeń z co najmniej trzema leptonami.');
  const cols = { 2: ['2', '(1,1), (2,0), (0,2)', '1', '1'], 3: ['3', '(2,1), (1,2)', '2', '1'], 4: ['4', '(2,2), (3,1), (1,3)', '4 lub 3', '1 + 1 lub 3'], 5: ['5+', 'np. (3,2)', '6+', '3 + 1 i więcej'] };
  const rows = [['n', 'Układ (e⁺, e⁻)', { t: 'Hipotezy\ne⁺e⁻', ...num }, { t: 'Hipotezy\njednoim.', ...num }, { t: 'Udział\nw tle [%]', ...num }, { t: 'Zostaje\ntła e⁺e⁻', ...num }, { t: 'Zostaje\ne⁺e⁺', ...num }, { t: 'Zostaje\ne⁻e⁻', ...num }]];
  T.classes.forEach(r => { const c = cols[r.c]; rows.push([c[0], c[1], { t: c[2], ...num }, { t: c[3], ...num }, { t: r.share, ...num }, { t: r.keptB, ...num, bold: true }, { t: r.keptP, ...num }, { t: r.keptM, ...num }]); });
  tbl(s, rows, { x: MX, y: 1.55, w: 8.5, colW: [0.55, 1.95, 0.95, 1.2, 1.0, 0.95, 0.95, 0.95], rowH: [0.62, 0.42, 0.42, 0.42, 0.42], fontSize: 12.5, margin: [2, 4, 2, 4] });
  cap(s, `n: liczba leptonów. M_{ee} od 0,14 do 0,80 ${M}, bez cięcia OA. Klasa 5+ ma tylko 66 par tła.`, MX, 4.12, 8.5, 0.3);
  fig(s, 'figs/classes.png', 9.3, 1.5, 3.43, 2.75, 'Ułamek par, które zostają, w funkcji liczby leptonów w zdarzeniu: tło e⁺e⁻ spada jak jeden przez liczbę hipotez, pary e⁺e⁺ i e⁻e⁻ zostają prawie wszystkie aż do czterech leptonów.', { left: true });
  cap(s, 'Rysunek: ile par zostaje po isBest == 1, zależnie od liczby leptonów w zdarzeniu.', 9.3, 3.8, 3.43, 0.5);
  bl(s, [
    'Klasa n = 2 nie traci nic i daje 56 % tła: to zdarzenia bez rywalizacji.',
    'Wzrost stosunku wzoru do tła pochodzi wyłącznie ze zdarzeń z trzema i więcej leptonami, czyli z 44 % tła.',
    'Pary jednoimienne zostają prawie wszystkie, bo w tych układach są zwykle pojedyncze.',
  ], { x: MX, y: 4.48, w: CW, h: 1.55 });
  result(s, 'tło e⁺e⁻ znika tylko w zdarzeniach z wieloma leptonami, a wzór tego nie uwzględnia.');
}

{
  const s = cs(P3('Trafność wyboru'), 'W zdarzeniach z rywalizacją isBest nie wyróżnia pary prawdziwej: zostaje jej tyle samo co tła',
    'To jest bezpośredni test tego, czy isBest wybiera parę prawdziwą. Gdyby wybierał ją trafnie, para prawdziwa zostawałaby częściej niż para tła. Tymczasem w zdarzeniach z trzema leptonami zostaje jej 0,515 w przedziale 0,14–0,80 GeV/c², tyle samo co tła, a poniżej 0,14 GeV/c² 0,494 wobec 0,496 dla tła. W zdarzeniach z czterema leptonami oba ułamki są bliskie 1/4. Wysokie 99,6 % z poprzednich slajdów bierze się stąd, że 99 % par prawdziwych leży w zdarzeniach z dwoma leptonami, gdzie nie ma z czego wybierać. χ² czasów przelotu nie rozróżnia tych par: wartości poniżej 1 ma 84–86 % par prawdziwych i 87–90 % par tła złożonego z prawdziwych leptonów. Powód leży w samym kryterium i pokazuje go następny slajd.');
  tbl(s, [
    ['Klasa', 'Przedział [GeV/c^{2}]', { t: 'Zostaje\ntła e⁺e⁻', ...num }, { t: 'Zostaje\npar prawdziwych', ...num }],
    ['n = 3', '0,00–0,14', { t: '0,496 ± 0,003 (37 768)', ...num }, { t: '0,494 ± 0,002 (71 534)', ...num }],
    ['n = 3', '0,14–0,80', { t: '0,515 ± 0,003 (37 417)', ...num }, { t: '0,515 ± 0,010 (2 322)', ...num }],
    ['n = 4', '0,00–0,14', { t: '0,249 ± 0,007 (3 633)', ...num }, { t: '0,261 ± 0,005 (8 487)', ...num }],
    ['n = 4', '0,14–0,80', { t: '0,246 ± 0,006 (5 576)', ...num }, { t: '0,208 ± 0,029 (202)', ...num }],
  ], { x: MX, y: 1.55, w: 7.1, colW: [0.75, 1.75, 2.3, 2.3], rowH: 0.38, fontSize: 13, ruleRows: [2] });
  cap(s, 'W nawiasach liczba par przed wyborem. Pary prawdziwe łącznie z parami „jedna cząstka dwa razy”. Bez cięcia na kąt otwarcia, niepewności statystyczne.', MX, 3.92, 7.1, 0.5);
  bl(s, [
    'Wysokie 99,6 % bierze się stąd, że 99 % par prawdziwych leży w zdarzeniach z dwoma leptonami, gdzie nie ma wyboru.',
    'Przyczynę widać w definicji kryterium w kodzie PAT3: następny slajd.',
  ], { x: MX, y: 4.48, w: 7.1, h: 1.55 });
  stat(s, 'około 1 / N', 'tyle par zostaje, zarówno prawdziwych, jak i tła; N to liczba hipotez e⁺e⁻ w zdarzeniu', 8.05, 1.55, 4.68, 1.3);
  stat(s, '84–86 % i 87–90 %', 'par prawdziwych i par tła złożonego z leptonów ma χ² czasów poniżej 1 (masa poniżej 0,14 GeV/c^{2}, OA > 2°)', 8.05, 3.0, 4.68, 1.3);
  tx(s, 'Dotyczy symulacji i obecnego kryterium wyboru. W danych nie ma historii śladów, ale jest ta sama preferencja detektora (następny slajd).', { x: 8.05, y: 4.5, w: 4.68, h: 1.4 });
  result(s, 'gdy w zdarzeniu jest kilka par, para prawdziwa nie zostaje częściej niż para tła.');
}

{
  const s = cs(P3('Powód'), 'Dlaczego isBest nie wskazuje pary prawdziwej: kryterium χ² porównuje tylko czasy przelotu',
    'Kod PAT3 (git.gsi.de/przygoda/pat3, master, commit 1cf58ad): HTimeCut::calcTof w user/htimecut.cc liczy χ² hipotezy, a HPidPool::fill w core/hpidpool.cc w każdym drzewie osobno daje isBest == 1 hipotezie o najmniejszym χ². Dla pary rachunek upraszcza się do wzoru na slajdzie. D to różnica zmierzonych czasów przelotu minus różnica czasów oczekiwanych z pędu, długości toru i masy elektronu. Wzór sprawdzono na wszystkich 8 039 114 parach tej symulacji: odtwarza zapisane χ² co do 0,1 %. Wzór nie zawiera żadnej informacji o tym, czy dwa leptony pochodzą z jednego rozpadu. Mówi tylko, czy oba są zgodne z prędkością światła i z jednym czasem startu. Dwa leptony z jednego zderzenia spełniają to tak samo, niezależnie od pochodzenia: D ma w zdarzeniach z trzema leptonami szerokość 0,167 ns dla par prawdziwych i 0,164 ns dla tła z leptonów. Rozróżnia tylko fałszywe ślady, dla których D jest przesunięte i szersze, ale te zostały już w większości usunięte przez cięcie PID. Wybór między hipotezami decyduje więc fluktuacja pomiaru czasu oraz detektor: w kodzie nominalna rozdzielczość RPC to 0,25 ns, a TOF 0,15 ns, więc hipoteza ze śladem w RPC ma mniejsze χ² i wygrywa w 70 % przypadków. W zdarzeniach, w których oba niewspólne ślady są w tym samym detektorze, para prawdziwa wygrywa w 51 %, czyli tyle samo co para tła. Gdy para prawdziwa ma ślad w RPC, wygrywa w 71 %, a gdy jej rywal ma ślad w RPC, tylko w 30 %. Ta sama preferencja jest w danych: przy kryterium kąta θ < 45° hipoteza z RPC wygrywa w 63,1 % w danych i w 63,0 % w symulacji. Bezpośredni test w danych: w zdarzeniach z trzema leptonami para o masie poniżej 0,03 GeV/c², czyli w dużej większości para prawdziwa, wygrywa z rywalem z przedziału 0,14–0,80 GeV/c² w 50,8 % z 462 223 zdarzeń. Rzeczywista rozdzielczość RPC jest lepsza od TOF, więc te wartości nominalne są odwrócone względem pomiaru; nie zmienia to jednak głównego wniosku, bo nawet przy poprawnych rozdzielczościach χ² nie mówi nic o wspólnym pochodzeniu pary.');
  const lw = 5.75;
  eq(s, 'chi2', MX, 1.6, 1.25, 'Chi kwadrat równa się D kwadrat przez cztery razy suma odwrotności kwadratów rozdzielczości, a D to różnica czasów zmierzonych minus różnica czasów oczekiwanych z pędu.');
  bl(s, [
    '**Prostymi słowami:** χ² sprawdza tylko, czy czasy przelotu dwóch leptonów pasują do jednego zderzenia. Nie wie, czy leptony pochodzą z jednego rozpadu.',
    '**Symbole:** t_{i} czas zmierzony, t_{i}^{p} oczekiwany z pędu; σ = 0,25 ns (RPC), 0,15 ns (TOF).',
    '**isBest == 1:** dostaje hipoteza o najmniejszym χ².',
    `**Bez informacji o pochodzeniu:** szerokość D to ${T.D.gen} ns dla par prawdziwych i ${T.D.lep} ns dla tła z leptonów.`,
    `**Decyduje detektor:** hipoteza ze śladem w RPC wygrywa w ${pct(T.choice_rpc.sim)}, bo większe σ daje mniejsze χ².`,
    `**To samo w danych:** para o masie poniżej 0,03 ${M} wygrywa z parą z 0,14–0,80 w ${pct(T.tag.all_d)} z ${T.tag.n} zdarzeń, a ślad w RPC w ${pct(T.choice_rpc.data)} (symulacja ${pct(T.choice_rpc.simT)}).`,
  ], { x: MX, y: 2.42, w: lw, h: 3.62, fontSize: 13, paraSpaceAfter: 4 });
  fig(s, 'figs/choice.png', MX + lw + 0.2, 1.55, CW - lw - 0.2, 4.4, 'Po lewej rozkład D dla par prawdziwych, tła z leptonów i tła z fałszywym śladem, zdarzenia z trzema leptonami: pierwsze dwa pokrywają się. Po prawej prawdopodobieństwo wyboru pary prawdziwej: 0,51 gdy oba niewspólne ślady są w tym samym detektorze, 0,71 gdy prawdziwa ma ślad w RPC, 0,30 gdy rywal.');
  result(s, 'χ² mówi, czy ślady są zgodne z jednym zderzeniem, a nie, czy pochodzą z jednego rozpadu.');
}

{
  const s = cs(P3('Błąd wzoru'), 'Skutek: przy isBest == 1 wzór przeszacowuje tło o 21 %',
    'Czarne punkty to stosunek wzoru do tła prawdziwego dla wszystkich kombinacji, czerwone dla isBest == 1. Różnica między nimi rozkłada się na dwa czynniki. Pierwszy, f jedn, mówi, jaka część wartości wzoru zostaje po isBest == 1: około 0,995. Drugi, f CB, mówi, jaka część tła prawdziwego e⁺e⁻ zostaje po isBest == 1: około 0,77. Stosunek wzoru do tła rośnie więc o czynnik 0,995 / 0,773. Równanie odtwarza wartość dla isBest == 1 bez dodatkowych założeń; jest to tożsamość, która pokazuje, skąd bierze się różnica. Wzór nie wie, że tło e⁺e⁻ zostało przycięte, bo przycięcie nie dotyka par jednoimiennych. Niezależnie od tego także przy isBest ≥ 0 wzór jest za mały poniżej 0,30 GeV/c²: R wynosi 0,908 w 0,14–0,30 GeV/c². To złamanie warunku akceptancji omówione na slajdzie 15; przy isBest == 1 ten sam błąd jest zamaskowany przez przeszacowanie, a stosunek 1,180 jest mniejszy niż 1,254 w 0,30–0,56 GeV/c².');
  eq(s, 'rdec', MX, 1.58, 2.2, 'R przy isBest == 1 równa się R przy isBest większym lub równym 0 razy f jedn podzielone przez f CB.');
  const rows = [['Przedział [GeV/c^{2}]', { t: 'R, isBest ≥ 0', ...num }, { t: 'f_{jedn}', ...num }, { t: 'f_{CB}', ...num }, { t: 'R, isBest == 1', ...num }]];
  T.dec.forEach(r => rows.push([r.w, { t: r.RA, ...num }, { t: r.fL, ...num }, { t: r.fB, ...num }, { t: r.RB, ...num, bold: true }]));
  tbl(s, rows, { x: MX, y: 2.45, w: 6.2, colW: [1.9, 1.2, 0.95, 0.95, 1.2], rowH: 0.34, fontSize: 14, ruleRows: [3], boldRows: [4] });
  bl(s, [
    '**R:** tło ze wzoru / tło prawdziwe.',
    '**f_{jedn}:** część wartości wzoru po isBest == 1.',
    '**f_{CB}:** część tła prawdziwego e⁺e⁻ po isBest == 1.',
    `**Także przy isBest ≥ 0** w 0,14–0,30 ${M} R = 0,908: złamany warunek akceptancji (slajd 15).`,
  ], { x: MX, y: 4.62, w: 6.2, h: 1.42, fontSize: 14, paraSpaceAfter: 3 });
  fig(s, 'figs/ratio_best.png', 7.05, 1.55, 5.68, 4.4, 'Stosunek tła ze wzoru do tła prawdziwego w funkcji masy pary dla isBest większego lub równego zero (czarne, blisko 1) i isBest == 1 (czerwone, od 1,1 do 1,6).');
  cap(s, 'Rysunek: ile razy tło ze wzoru różni się od tła prawdziwego w obu wariantach; 1 oznacza zgodność.', 7.05, 5.55, 5.68, 0.4);
  result(s, '0,938 × 0,995 / 0,773 = 1,207: wzór prawie się nie zmienia, a tło traci 23 %.');
}

{
  const s = cs(P3('Skutek dla sygnału'), 'Przy isBest == 1 sygnał jest zaniżony, a błąd rośnie wraz z masą',
    'Sygnał ze wzoru to wszystkie pary minus tło ze wzoru. Przy isBest ≥ 0 stosunek do sygnału prawdziwego leży blisko jedności. Przy isBest == 1 sygnał jest zaniżony tym bardziej, im większa masa: 2,2 % w przedziale 0,05–0,14 GeV/c², 3,7 % w 0,14–0,30 i 12,4 % w 0,30–0,56. Błąd rośnie, bo stosunek sygnału do tła maleje z masą, a nadmiar tła ze wzoru jest odejmowany od coraz mniejszego sygnału.');
  const rows = [['Przedział [GeV/c^{2}]', { t: 'isBest ≥ 0', ...num }, { t: 'isBest == 1', ...num }]];
  T.windows.filter(r => r.w !== '0,00–0,02' && r.w !== '0,02–0,05').forEach(r => rows.push([r.w, { t: r.all.Sr, ...num }, { t: r.best.Sr, ...num, bold: true }]));
  tbl(s, rows, { x: MX, y: 1.55, w: 5.2, colW: [2.2, 1.5, 1.5], rowH: 0.4, fontSize: 15, ruleRows: [3], boldRows: [4] });
  cap(s, 'Sygnał ze wzoru (wszystkie − wzór) podzielony przez sygnał prawdziwy, każdy wariant względem własnego sygnału prawdziwego.', MX, 3.65, 5.2, 0.7);
  bl(s, [
    '**Rysunek:** u góry sygnał, na dole sygnał ze wzoru podzielony przez prawdziwy.',
    'Błąd rośnie wraz z masą, bo stosunek sygnału do tła maleje.',
  ], { x: MX, y: 4.35, w: 5.2, h: 1.6 });
  fig(s, 'figs/signal_best.png', 6.15, 1.5, 6.58, 4.5, 'Sygnał w funkcji masy pary: sygnał prawdziwy, wszystkie pary minus wzór dla isBest większego lub równego 0 i dla isBest == 1. Dolny panel: stosunek do sygnału prawdziwego.');
  result(s, `sygnał przy isBest == 1 jest zaniżony o 5,5 % w 0,14–0,80 ${M} i o 12 % w 0,30–0,56 ${M}.`);
}

{
  const s = cs(P3('Wszystkie kombinacje'), 'Do czego prowadzi isBest ≥ 0: wzór zgadza się z tłem, a istotność spada tylko o 2 %',
    'Oba panele mają tę samą skalę. Pary e⁺e⁺ i e⁻e⁻ są w obu praktycznie takie same, dlatego wzór, czerwone kwadraty, też się nie zmienia. Zmienia się tło prawdziwe e⁺e⁻, czarna linia: przy isBest ≥ 0 leży na kwadratach, a przy isBest == 1 pod nimi. Tabela zbiera to dla przedziału 0,14–0,80 GeV/c². Wybór wszystkich kombinacji dodaje 23 687 par e⁺e⁻, z czego 22 400 to tło, a tylko 1287 to sygnał. Te dodatkowe pary tła są jednak opisane przez wzór, bo w parach jednoimiennych liczba kombinacji jest pełna: stosunek wzoru do tła prawdziwego wynosi 0,938 ± 0,004 zamiast 1,207 ± 0,006, a sygnał ze wzoru jest poprawny do 2,1 % zamiast 5,5 %. Dolne paski pokazują stosunek tła prawdziwego do tła ze wzoru, czyli wprost poprawkę k: przy isBest ≥ 0 powyżej 1 w 0,14–0,30 GeV/c², a przy isBest == 1 około 0,8 w całym zakresie. Pozostała różnica 6 % tła to poprawka akceptancji k, wyznaczona w części czwartej. Kosztem jest niewielka utrata istotności statystycznej, czyli S/√(N + σ²): z 425 do 416, o 2 %.');
  fig(s, 'figs/cb_both.png', MX, 1.5, 7.35, 4.5, 'Dwa panele o tej samej skali: isBest większe lub równe 0 oraz isBest == 1. W każdym tło prawdziwe e⁺e⁻, tło ze wzoru, pary e⁺e⁺ i e⁻e⁻ w funkcji masy pary. Pod spodem stosunek tła prawdziwego do tła ze wzoru: przy isBest ≥ 0 blisko 1, poza zakresem 0,14–0,30 GeV/c²; przy isBest == 1 około 0,8.', { left: true });
  const a = T.compare.all, b = T.compare.best, row = (t, k) => [t, { t: a[k], ...num }, { t: b[k], ...num }];
  tbl(s, [
    ['0,14–0,80 GeV/c^{2}', { t: 'isBest ≥ 0', ...num }, { t: 'isBest == 1', ...num }],
    row('Pary e⁺e⁻', 'N'), row('Sygnał prawdziwy', 'S'), row('Tło prawdziwe', 'B'), row('Pary e⁺e⁺', 'pp'), row('Pary e⁻e⁻', 'mm'), row('Tło ze wzoru', 'CB'),
    ['Wzór / tło prawdziwe', { t: `${a.R} ± ${a.eR}`, ...num, bold: true }, { t: `${b.R} ± ${b.eR}`, ...num, bold: true }],
    ['Sygnał ze wzoru / prawdziwy', { t: a.Sr, ...num, bold: true }, { t: b.Sr, ...num, bold: true }],
    row('Istotność', 'fom'),
  ], { x: 8.05, y: 1.5, w: 4.68, colW: [2.1, 1.29, 1.29], rowH: 0.38, fontSize: 13, ruleRows: [6, 8] });
  cap(s, 'Dolne paski: tło prawdziwe podzielone przez tło ze wzoru. Istotność: S/√(N + σ^{2}).', 8.05, 5.55, 4.68, 0.45);
  result(s, 'przy isBest ≥ 0 sygnał jest poprawny do 2,1 %, a koszt to tylko 2 % istotności (425 → 416).');
}

// ================================================================ część 4
divider('Część 4: poprawka k', 'Część 4: poprawka k z mieszania zdarzeń w danych',
  ['Procedura: mieszanie zdarzeń w danych', 'Test metody na symulacji i znaczenie składu puli', 'k w danych i sprawdzenia'],
  'W czwartej części wyznaczono poprawkę k, która w poprzednich częściach była przyjęta jako równa 1. Poprawka pochodzi wyłącznie z danych. Symulacja służy tylko do sprawdzenia samej metody mieszania, bo tła kombinatorycznego w danych nie odtwarza.');

{
  const s = cs(P4('Procedura'), 'Jak wyznaczyć poprawkę k: mieszanie zdarzeń w danych',
    'Procedura ma siedem kroków i korzysta wyłącznie z danych. Dotyczy wszystkich kombinacji, czyli isBest ≥ 0. Krok 1: zdarzenia z trygerem PT3, po cięciu na wierzchołek i PID 5σ. Krok 2: z wierszy par e⁺e⁻ odtwarzamy niepowtarzalne ślady w każdym zdarzeniu. Krok 3: każde zdarzenie mieszamy z 10 poprzednimi zdarzeniami tej samej klasy: n = 2 albo n ≥ 3. Krok 4: k_mix liczymy osobno w każdym przedziale masy i każdej klasie. W mieszanych zdarzeniach nie ma par skorelowanych, więc k_mix mierzy różnicę akceptancji. Krok 5: wynikiem jest suma po klasach ważona ich udziałem we wzorze w danych. Krok 6: do punktów dopasowujemy gładką funkcję: rozkład Gaussa o różnej szerokości po obu stronach maksimum plus stała, a pasmo niepewności, ze statystyki i różnicy puli, przenosimy na sygnał. Krok 7: tło to k razy wzór. Poniżej 0,05 GeV/c² k nie opisuje akceptancji, bo dominują tam fałszywe ślady i jedna cząstka zrekonstruowana dwa razy; tam potrzebne jest cięcie na kąt otwarcia. Symulacja nie wchodzi do k: nie odtwarza tła kombinatorycznego w danych. Skrypt: plots/k_mix_pools.C.');
  tx(s, paras([
    '**1. Te same zdarzenia:** wierzchołek, PID 5σ, wszystkie kombinacje (isBest ≥ 0), PT3.',
    '**2. Leptony:** z par e⁺e⁻ odtwarzamy niepowtarzalne ślady e⁺ i e⁻ w każdym zdarzeniu.',
    '**3. Mieszanie:** każde zdarzenie z 10 poprzednimi tej samej klasy: n = 2 albo n ≥ 3.',
    '**4. k_{mix,n}(M):** par skorelowanych nie ma, więc k_{mix} mierzy różnicę akceptancji.',
    '**5. Wynik:** suma po klasach ważona ich udziałem F_{n} we wzorze w danych.',
    '**6. Gładka funkcja:** dopasowanie do punktów k(M), pasmo niepewności przeniesione na sygnał.',
    '**7. Użycie:** tło = k(M) · 2√(N_{++} N_{−−}) przy isBest ≥ 0.',
  ], { fontSize: 15 }), { x: MX, y: 1.55, w: 7.2, h: 4.4, paraSpaceAfter: 7, name: 'Kroki' });
  card(s, 'Równania', '', 8.05, 1.55, 4.68, 4.4, { fill: C.accent1 });
  eq(s, 'kmix', 8.25, 2.1, 1.32, 'k mix równa się liczba par różnoimiennych w mieszaniu przez dwa pierwiastki z iloczynu par jednoimiennych w mieszaniu.');
  eq(s, 'kfinal', 8.25, 3.15, 1.3, 'k równa się jeden przez F razy suma po klasach n z F n razy k mix n, gdzie F to dwa pierwiastki z iloczynu liczb par jednoimiennych.');
  tx(s, `**Symulacja nie wchodzi do k:** nie odtwarza tła kombinatorycznego w danych. Poniżej 0,05 ${M} k nie opisuje akceptancji.`, { x: 8.25, y: 4.1, w: 4.3, h: 1.7, fontSize: 14, color: C.text2 });
  result(s, 'k pochodzi wyłącznie z danych, z mieszania zdarzeń w dwóch klasach krotności leptonów.');
}

{
  const s = cs(P4('k w symulacji'), `Test metody na symulacji: mieszanie odtwarza kształt k, ale nie pełną wartość`,
    'Czarne punkty to k_true. Liczymy je tylko w symulacji: w każdym przedziale masy tło prawdziwe policzone z historii śladów (zdarzenia PT3, isBest ≥ 0) dzielimy przez wzór 2√(N₊₊N₋₋) z tych samych zdarzeń. Wzór pomnożony przez k_true daje dokładnie tło prawdziwe tej próbki. Nie jest to k dla danych, tylko punkt odniesienia, z którym porównujemy mieszanie zdarzeń. Maksimum przypada w 0,14–0,20 GeV/c² i wynosi 1,144. Mieszanie zdarzeń odtwarza kształt, ale jest niższe: 1,060 dla puli wszystkich zdarzeń i 1,083 dla puli zdarzeń z co najmniej trzema leptonami. Różnica między pulami pokazuje, że wynik mieszania zależy od tego, jakie leptony trafią do puli. Pula ze wszystkich zdarzeń składa się głównie z leptonów par prawdziwych. Pula z n ≥ 3 jest bliższa leptonom tworzącym tło i bliższa wartości prawdziwej. Poniżej 0,05 GeV/c² poprawka nie opisuje akceptancji. Liczby dla zdarzeń PT3, tak jak w danych. Symulacja nie odtwarza tła kombinatorycznego w danych, dlatego to tylko test metody na próbce ze znaną prawdą, a różnicy między mieszaniem a k_true nie przenosimy na dane.');
  fig(s, 'figs/k_sim.png', MX, 1.5, 7.0, 4.5, 'k w funkcji masy pary w symulacji: k true z historii śladów oraz k z mieszania zdarzeń dla puli wszystkich zdarzeń i puli zdarzeń z co najmniej trzema leptonami. Maksimum około 1,15 przy 0,15 GeV/c².', { left: true });
  const rows = [['Przedział\n[GeV/c²]', { t: 'k_{true}', ...num }, { t: 'k_{mix}\nwszystkie', ...num }, { t: 'k_{mix}\nn ≥ 3', ...num }]];
  T.k.forEach(r => rows.push([r.w, { t: `${r.kt} ± ${r.ekt}`, ...num, bold: true }, { t: r.ksa, ...num }, { t: r.ksn, ...num }]));
  tbl(s, rows, { x: 7.85, y: 1.5, w: 4.88, colW: [1.13, 1.45, 1.2, 1.1], rowH: [0.56, 0.3, 0.3, 0.3, 0.3, 0.3], fontSize: 13 });
  card(s, 'k_{true}: punkt odniesienia, tylko w symulacji',
    'Tło prawdziwe ze śladów (część 2) podzielone przez wzór 2√(N_{++}N_{−−}), zdarzenia PT3, isBest ≥ 0. Wzór pomnożony przez k_{true} daje tło prawdziwe tej próbki. To nie jest k dla danych.',
    7.85, 3.66, 4.88, 1.58, { fill: C.accent1, headSize: 14, size: 13.5, gap: 2 });
  bl(s, [
    '**k_{mix}:** ten sam stosunek dla par z mieszania zdarzeń, czyli samego tła.',
    'Test metody: różnicy nie przenosimy na dane.',
  ], { x: 7.85, y: 5.28, w: 4.88, h: 0.78, fontSize: 13.5, paraSpaceAfter: 2 });
  result(s, `w 0,14–0,20 ${M} mieszanie daje w symulacji k_{mix} = 1,06–1,08 wobec k_{true} = 1,14.`);
}

{
  const s = cs(P4('Skład puli'), `Skąd różnica z k_{true}: wynik mieszania zależy od tego, jakie leptony są w puli`,
    'Mieszanie zdarzeń łączy leptony z różnych zdarzeń, a więc nie zna par skorelowanych, ale jego wynik zależy od tego, jakie leptony trafiły do puli. W puli ze wszystkich zdarzeń 95 % stanowią zdarzenia z jedną parą prawdziwą: oba leptony tego samego rozpadu zostały zarejestrowane, więc mają inne pędy i kąty niż leptony par tła, których partner zwykle przepadł. Dlatego taka pula daje najmniejsze k_mix. W tabeli pokazano przedział 0,14–0,30 GeV/c²; wiersze oznaczone jako wybór z historii śladów wymagają prawdy z symulacji, więc w danych są niedostępne. Prawdziwa wartość k_true w symulacji, 1,101, leży wewnątrz zakresu, ale żadna pula dostępna w danych jej nie odtwarza. Tej różnicy nie przenosimy na dane, bo symulacja nie odtwarza tła kombinatorycznego w danych. Pula nie tworzy osobnych poprawek: wynikiem jest jedna funkcja k(M). Powyżej 0,30 GeV/c² różnica znika. Makro: diag/k_pool_composition_sim.cc.');
  tbl(s, [
    ['Pula leptonów (symulacja, 0,14–0,30)', { t: 'k_{mix}', ...num }],
    ['leptony par prawdziwych, n = 2 (historia śladów)', { t: '1,031', ...num }],
    ['wszystkie zdarzenia', { t: '1,035', ...num }],
    ['zdarzenia z n ≥ 3', { t: '1,056', ...num }],
    ['n = 2 z parą tła (historia śladów)', { t: '1,100', ...num }],
    ['n ≥ 3 oraz zdarzenia tylko z parami jednoimiennymi', { t: '1,137', ...num }],
    ['k_{true} (tło prawdziwe / wzór)', { t: '1,101', ...num, bold: true }],
  ], { x: MX, y: 1.55, w: 7.1, colW: [5.9, 1.2], rowH: 0.48, fontSize: 15, boldRows: [6], ruleRows: [5] });
  cap(s, 'Zdarzenia z n ≥ 3 i wszystkie zdarzenia to pule możliwe do zbudowania w danych. Pozostałe są testem w symulacji.', MX, 5.15, 7.1, 0.5);
  card(s, 'Co z tego wynika', [
    '**Jedno k:** pule nie dają dwóch poprawek. n = 2 i n ≥ 3 to dwie klasy, ważone udziałem w danych.',
    '**Różnica z k_{true}:** pochodzi głównie ze składu puli i znika powyżej 0,30 GeV/c^{2}.',
    '**Granica metody:** żadna pula nie odtwarza k_{true} w symulacji. Tej różnicy nie przenosimy na dane.',
  ], 8.05, 1.55, 4.68, 4.4, { size: 15, gap: 10, fill: C.accent1 });
  result(s, 'k_{mix} zmienia się od 1,03 do 1,14 zależnie od puli, a k_{true} = 1,10 leży wewnątrz tego zakresu.');
}

{
  const s = cs(P4('k w danych'), `Dane: k = ${T.k_int.k} ± ${T.k_int.err} w 0,14–0,80 ${M}, a maksimum to ${T.k[1].kw} w 0,14–0,20 ${M}`,
    `Czarne punkty to k z mieszania zdarzeń w danych w każdym przedziale masy: klasy n = 2 i n ≥ 3 ważone ich udziałem we wzorze w danych. Czerwona krzywa to gładkie dopasowanie: rozkład Gaussa o różnej szerokości po obu stronach maksimum i stała (χ² = ${T.kfit.chi2} na ${T.kfit.ndf} stopni swobody, szerokość lewego zbocza ustalona). Pasmo obejmuje statystykę oraz różnicę między pulami. Wartość w przedziale 0,14–0,80 GeV/c² to średnia ważona tłem w danych. Sygnał w 0,14–0,80 GeV/c² wynosi ${D4.Sk} ± ${D4.eS} par, a bez poprawki ${D4.Sa}. Symulacja nie wchodzi do tej poprawki. Gdyby jednak użyć stosunku k_true do k_mix z symulacji, sygnał wyniósłby ${D4.Skm} par; nie stosujemy tego, bo symulacja nie odtwarza tła w danych, ale jest to miara, jak bardzo wynik zależałby od takiego założenia.`);
  fig(s, 'figs/k_data.png', MX, 1.5, 7.0, 4.5, 'k w funkcji masy pary w danych: punkty z klas krotności, gładkie dopasowanie z pasmem niepewności oraz k z mieszania dla dwóch pul.', { left: true });
  const rows = [['Przedział\n[GeV/c²]', { t: 'k_{mix}^{dane}\nwszystkie', ...num }, { t: 'k_{mix}^{dane}\nn ≥ 3', ...num }, { t: 'k', ...num }]];
  T.k.forEach(r => rows.push([r.w, { t: r.kda, ...num }, { t: r.kdn, ...num }, { t: `${r.kw} ± ${r.ekw}`, ...num, bold: true }]));
  tbl(s, rows, { x: 7.85, y: 1.55, w: 4.88, colW: [1.13, 1.2, 1.05, 1.5], rowH: [0.62, 0.38, 0.38, 0.38, 0.38, 0.38], fontSize: 13.5 });
  bl(s, [
    '**Rysunek:** k z mieszania zdarzeń w danych; czerwona krzywa to wygładzenie, którego używamy.',
    'Poprawka pochodzi z samych danych; pasmo obejmuje statystykę i różnicę między pulami.',
  ], { x: 7.85, y: 4.38, w: 4.88, h: 1.62, fontSize: 14 });
  result(s, `k zmienia sygnał w 0,14–0,80 o ${num2(D4.SkN/D4.SaN)} %: ${D4.Sa} → ${D4.Sk} ± ${D4.eS} par.`);
}

{
  const s = cs(P4('Sprawdzenia'), 'Sprawdzenia k w danych: brak zależności od trygera i mała zależność od pędu poprzecznego',
    'Dwa sprawdzenia, oba na danych. Tryger: k z mieszania policzone osobno dla PT2 i PT3. Poniżej 0,56 GeV/c² różnica wynosi najwyżej 0,4 %, więc PT3 można zostawić: efekt trygera zmienia liczbę par, ale w stosunku k się skraca. Pozostaje on w poprawce sygnału na tryger. Pęd poprzeczny pary: k z mieszania liczone w siatce masa i p_T, ważone rozkładem w danych, porównano z wersją bez p_T. Zmiana wynosi co najwyżej 2 % przy 0,14–0,16 GeV/c² i 0,3 % powyżej 0,20 GeV/c².');
  const w1 = 5.6, w3 = 5.6, g_ = CW - w1 - w3;
  tx(s, '**Tryger: k z mieszania w danych**', { x: MX, y: 1.55, w: w1, h: 0.3, fontSize: 16 });
  const r1 = [['Przedział [GeV/c²]', { t: 'PT2', ...num }, { t: 'PT3', ...num }]]; T.trig.forEach(r => r1.push([r.w, { t: r.k2, ...num }, { t: r.k3, ...num, bold: true }]));
  tbl(s, r1, { x: MX, y: 1.9, w: w1, colW: [2.6, 1.5, 1.5], rowH: 0.4, fontSize: 15 });
  tx(s, '**p_{T}: k z siatki masa–p_{T} / k bez p_{T}**', { x: MX + w1 + g_, y: 1.55, w: w3, h: 0.3, fontSize: 16 });
  const r3 = [['Przedział [GeV/c²]', { t: 'pula wszystkich', ...num }, { t: 'pula n ≥ 3', ...num }]]; T.cpt.forEach(r => r3.push([r.w, { t: r.all, ...num }, { t: r.n3, ...num }]));
  tbl(s, r3, { x: MX + w1 + g_, y: 1.9, w: w3, colW: [2.2, 1.8, 1.6], rowH: 0.34, fontSize: 14 });
  bl(s, [
    `**Tryger:** PT3 zostaje. Poniżej 0,56 ${M} różnica z PT2 to najwyżej 0,4 %; efekt trygera na sam sygnał jest osobną poprawką.`,
    `**Wagi ε_{n}:** w PT2 i PT3 różnią się najwyżej o 1,3 % dla n = 3 i o 3 % dla n = 4.`,
    `**p_{T}:** zmiana do 2 % przy 0,14–0,16 ${M}, powyżej 0,20 ${M} do 0,3 %, więc dla widma masy k(M) wystarcza.`,
  ], { x: MX, y: 4.72, w: CW, h: 1.32, fontSize: 14, paraSpaceAfter: 3 });
  result(s, 'k nie zależy od trygera, a zależność od p_{T} jest mała i ogranicza się do najmniejszych mas.');
}

// ================================================================ część 5
divider('Część 5: dane', 'Część 5: dane eksperymentalne',
  ['Widmo, tło ze wzoru i sygnał przy isBest ≥ 0 oraz isBest == 1', 'Ile par zostaje po isBest == 1', 'Skąd wagi ε_{n} dla isBest == 1', 'Sygnał w obu wariantach', 'k i ε_{n} w innych rozkładach'],
  'W piątej części te same wielkości pokazano dla danych eksperymentalnych: p+p, dni 034–067, zdarzenia PT3, po cięciu na wierzchołek i cięciach PID dla danych. Historii śladów tu nie ma.');

{
  const s = cs(P5('Widma'), 'Dane: widmo e⁺e⁻, tło ze wzoru i sygnał przy isBest ≥ 0 i przy isBest == 1',
    'Czarne punkty to wszystkie pary e⁺e⁻, czerwone tło ze wzoru, niebieskie sygnał. Lewy panel to wszystkie kombinacje, prawy jedna para na zdarzenie. Obie osie są takie same. Tło ze wzoru prawie się nie zmienia, a wszystkich par e⁺e⁻ ubywa, więc różnica, czyli sygnał, jest wyraźnie mniejsza po prawej stronie. Dolne paski pokazują stosunek sygnału do tła. W danych tło przeważa nad sygnałem w szerokim zakresie mas, od około 0,1 do 0,6 GeV/c². Przy 0,78 GeV/c² widać ω.');
  fig(s, 'figs/data_spectra.png', MX, 1.5, 8.7, 4.5, 'Dwa panele: isBest większe lub równe 0 oraz isBest == 1. Czarne punkty: wszystkie pary e⁺e⁻; czerwone: tło ze wzoru; niebieskie: sygnał. Pod spodem stosunek sygnału do tła.', { left: true });
  tbl(s, [
    ['Pary', { t: 'isBest ≥ 0', ...num }, { t: 'isBest == 1', ...num }],
    ['e⁺e⁻', { t: '47 035 765', ...num }, { t: '45 042 282', ...num }],
    ['e⁺e⁺', { t: '1 772 091', ...num }, { t: '1 712 925', ...num }],
    ['e⁻e⁻', { t: '5 348 663', ...num }, { t: '4 949 874', ...num }],
  ], { x: 9.5, y: 1.55, w: 3.23, colW: [0.83, 1.2, 1.2], rowH: 0.38, fontSize: 13 });
  bl(s, [
    '**Rysunek:** po lewej wszystkie kombinacje, po prawej jedna para na zdarzenie. Czarne to wszystkie pary, czerwone tło, niebieskie sygnał.',
    'Dane p+p, zdarzenia PT3, po cięciu na wierzchołek.',
    'Tło ze wzoru prawie się nie zmienia, a par e⁺e⁻ ubywa.',
    `Tło przeważa od 0,1 do 0,6 ${M}.`,
  ], { x: 9.5, y: 3.25, w: 3.23, h: 2.75, fontSize: 14 });
  result(s, 'tło przeważa nad sygnałem w szerokim zakresie mas, więc błąd wzoru silnie wpływa na sygnał.');
}

{
  const d = T.data[4];
  const s = cs(P5('Co zostaje'), `Dane: isBest == 1 usuwa ${(100*(1 - parseFloat(d.eps.replace(',', '.')))).toFixed(0)} % par e⁺e⁻, a tylko ${(100*(1 - parseFloat(d.rC.replace(',', '.')))).toFixed(0)} % tła ze wzoru`,
    'Wykres pokazuje stosunek liczby par z isBest == 1 do liczby wszystkich par, osobno dla e⁺e⁻, e⁺e⁺ i e⁻e⁻, oraz stosunek tła ze wzoru w obu wariantach. Ten sam mechanizm działa w danych i w symulacji: pary e⁺e⁻ tracą znacznie więcej niż pary jednoimienne, bo rywalizują w każdym zdarzeniu z wieloma leptonami. W przedziale 0,14–0,80 GeV/c² znika 17 % par e⁺e⁻, a tło ze wzoru maleje tylko o 5,3 %. Pary e⁻e⁻ tracą więcej niż e⁺e⁺, bo w danych jest ich trzykrotnie więcej. W przeciwieństwie do symulacji nie można tu rozdzielić sygnału i tła na podstawie historii śladów.');
  fig(s, 'figs/data_kept.png', MX, 1.5, 6.9, 4.5, 'Stosunek liczby par z isBest == 1 do wszystkich par w funkcji masy, osobno dla e⁺e⁻, e⁺e⁺ i e⁻e⁻, oraz stosunek tła ze wzoru. Pary e⁺e⁻ spadają do około 0,82, pary jednoimienne i wzór leżą około 0,93–0,97.', { left: true });
  const rows = [['Przedział\n[GeV/c²]', { t: 'e⁺e⁻', ...num }, { t: 'e⁺e⁺', ...num }, { t: 'e⁻e⁻', ...num }, { t: 'tło ze\nwzoru', ...num }]];
  T.data.forEach(r => rows.push([r.w, { t: r.eps, ...num, bold: true }, { t: r.epp, ...num }, { t: r.emm, ...num }, { t: r.rC, ...num, bold: true }]));
  tbl(s, rows, { x: 7.75, y: 1.55, w: 4.98, colW: [1.42, 0.89, 0.89, 0.89, 0.89], rowH: [0.62, 0.38, 0.38, 0.38, 0.38, 0.38], fontSize: 14, ruleRows: [4], boldRows: [5] });
  cap(s, 'Stosunek isBest == 1 do isBest ≥ 0. Dla tła ze wzoru: 2√(N_{++} N_{−−}) policzone z par w obu wariantach.', 7.75, 4.32, 4.98, 0.5);
  tx(s, '**Rysunek:** jaka część par zostaje po isBest == 1. Mechanizm jest taki sam jak w symulacji: pary e⁺e⁻ częściej mają rywala.', { x: 7.75, y: 4.92, w: 4.98, h: 1.05, fontSize: 15 });
  result(s, 'ta sama asymetria co w symulacji: pary e⁺e⁻ tracą kilka razy więcej niż tło ze wzoru.');
}

{
  const s = cs(P5('Wagi εₙ'), 'Skąd wagi klas n: podmieniamy lepton na lepton z innego zdarzenia i patrzymy, czy nowa para wygrywa',
    'Wagi ε_n pochodzą w całości z danych. Po co są: po wyborze isBest == 1 w zdarzeniu z n leptonami zostaje tylko część par tła, a zwykły wzór tego nie wie, bo pary jednoimienne zwykle nie mają rywali. Procedura ma cztery kroki. Bierzemy prawdziwe zdarzenie. Jeden z leptonów zastępujemy leptonem tego samego znaku z wcześniejszego zdarzenia, dzięki czemu nowe pary są na pewno parami tła, a reszta zdarzenia zostaje ta sama. Dla nowych par liczymy χ² i identyfikację tak samo jak PAT i układamy ranking razem z pozostałymi parami zdarzenia. Ułamek nowych par, które wygrały, to ε_n. Powtarzamy to dla każdego leptonu w zdarzeniu. Wartości są bliskie jeden przez liczbę hipotez e⁺e⁻. Tło przy isBest == 1 liczymy następnie w każdej klasie n osobno: ε_n razy dwa pierwiastki z iloczynu liczb par jednoimiennych z danej klasy, a sumę mnożymy przez k. Dla par jednoimiennych tę samą metodą można przewidzieć ułamek wygranych i porównać z ułamkiem zmierzonym: zgadzają się do 4 %. Makro: plots/isbest_emulation_betanew.C.');
  const steps = [['1. Bierzemy zdarzenie', 'Prawdziwe zdarzenie z danych, z n leptonami.'],
    ['2. Podmieniamy lepton', 'Jeden lepton zastępujemy leptonem tego samego znaku z innego zdarzenia. Nowe pary są na pewno parami tła.'],
    ['3. Układamy ranking', 'Dla nowych par liczymy χ² i identyfikację tak samo jak PAT i szeregujemy je z resztą par zdarzenia.'],
    ['4. Liczymy wygrane', 'Ułamek nowych par, które wygrały, to ε_{n}. Powtarzamy dla każdego leptonu.']];
  const sw = (CW - 3*0.15)/4;
  steps.forEach((st, i) => card(s, st[0], st[1], MX + i*(sw + 0.15), 1.55, sw, 1.8, { size: 14, headSize: 16, fill: C.accent1 }));
  const rows = [['Leptony\nn', { t: 'Hipotezy\ne⁺e⁻', ...num }, { t: 'ε_{n}\n0,05–0,14', ...num }, { t: 'ε_{n}\n0,14–0,30', ...num }, { t: 'ε_{n}\n0,30–0,80', ...num }, { t: 'Udział we\nwzorze [%]', ...num }]];
  T.eps.forEach(r => rows.push([r.n, { t: r.hyp, ...num }, { t: r.e[0], ...num }, { t: r.e[1], ...num, bold: true }, { t: r.e[2], ...num }, { t: r.share, ...num }]));
  tbl(s, rows, { x: MX, y: 3.5, w: 6.7, colW: [1.15, 1.1, 1.1, 1.1, 1.1, 1.15], rowH: [0.6, 0.36, 0.36, 0.36, 0.36], fontSize: 13 });
  const ls = T.eps_ls;
  bl(s, [
    '**Po co:** po isBest == 1 z tła zostaje tylko część par, a zwykły wzór tego nie wie.',
    `**Użycie:** tło w każdej klasie n mnożymy przez ε_{n}, sumujemy i mnożymy przez k. Tło spada do ${T.weightCB.r} wartości ze zwykłego wzoru.`,
    `**Sprawdzenie:** dla e⁺e⁺ przy n = 3 metoda daje ${ls.pp.n3.e}, pomiar ${ls.pp.n3.d}; dla e⁻e⁻ ${ls.mm.n3.e} i ${ls.mm.n3.d}.`,
  ], { x: 7.55, y: 3.5, w: 5.18, h: 2.5, fontSize: 14, paraSpaceAfter: 5 });
  result(s, 'ε_{n} jest bliskie 1 / (liczba hipotez e⁺e⁻) i pochodzi w całości z danych.');
}

{
  const d = T.data[4];
  const s = cs(P5('Sygnał'), `Dane, 0,14–0,80 ${M}: zwykły wzór przy isBest == 1 zaniża sygnał o ${num2(D4.SbN/D4.SkN)} %, z poprawnym tłem różnica to ${((1 - parseFloat(D4.re.replace(',', '.')))*100).toFixed(1).replace('.', ',')} %`,
    `Slajd porównuje trzy sposoby liczenia sygnału. Czarne punkty: S = N₊₋ − k · 2√(N₊₊N₋₋), wszystkie kombinacje, k to poprawka z części czwartej, szare pasmo to jej niepewność. Czerwone kwadraty: S = N₊₋ − 2√(N₊₊N₋₋), wszystko z par isBest == 1, czyli zwykły wzór. Zielone trójkąty: isBest == 1, ale tło liczone osobno w klasach krotności n: S = N₊₋ − k · Σₙ εₙ · 2√(N₊₊,ₙ N₋₋,ₙ). εₙ to ułamek par tła wygrywających rywalizację w klasie n, wyznaczony z emulacji rankingu na danych. Zwykły wzór zaniża sygnał w 0,14–0,80 GeV/c² o ${num2(D4.SbN/D4.SkN)} %. Korekta tła usuwa tylko tło, więc para prawdziwa, która przegrała ranking, nie wraca do próbki; dlatego zielone punkty leżą niżej niż czarne. Wariant z wagami klas zgadza się z czarnymi punktami w granicach ${((parseFloat(D4.re.replace(',', '.')) - 1)*100).toFixed(1).replace('.', ',')} %. Bez poprawki k wynik dla isBest ≥ 0 byłby o ${((D4.SaN/D4.SkN - 1)*100).toFixed(0)} % większy, więc k jest w danych równie ważne jak wybór wariantu. Wagi εₙ pochodzą z emulacji, a nie z prawdy. Rysunek ma przedziały 5 MeV/c² poniżej 0,14 GeV/c² i 20 MeV/c² powyżej; poniżej 0,05 GeV/c² przyjęto k = 1.`);
  fig(s, 'figs/data_signal.png', MX, 1.5, 6.4, 4.5, 'Sygnał w funkcji masy pary dla trzech wariantów z szarym pasmem niepewności k. Dolny panel: stosunek do wariantu isBest większe lub równe 0 z poprawką k. Zwykły wzór dla isBest == 1 leży poniżej, od 0,5 do 1; wariant z wagami klas przy jedności.', { left: true });
  const rows = [['Przedział\n[GeV/c²]', { t: 'isBest ≥ 0\nk = 1', ...num }, { t: 'isBest ≥ 0\nz k', ...num }, { t: 'isBest == 1\nzwykły wzór', ...num }, { t: 'isBest == 1\nwagi n, k', ...num }]];
  T.data.forEach(r => rows.push([r.w, { t: r.Sa, ...num }, { t: r.Sk, ...num, bold: true }, { t: r.Sb, ...num }, { t: r.Se, ...num, bold: true }]));
  tbl(s, rows, { x: 7.2, y: 1.55, w: 5.53, colW: [1.05, 1.06, 1.06, 1.18, 1.18], rowH: [0.62, 0.38, 0.38, 0.38, 0.38, 0.38], fontSize: 12, ruleRows: [4], boldRows: [5] });
  cap(s, `Liczba par sygnału. Niepewność z k w 0,14–0,80: ± ${D4.eS} par. Rysunek: przedziały 5 MeV/c^{2} do 0,14; poniżej 0,05 k = 1.`, 7.2, 4.3, 5.53, 0.45);
  tx(s, [...rt('**Zwykły wzór:** N_{+−} − 2√(N_{++} N_{−−}). **Wagi n, k:** N_{+−} − k · Σ_{n} ε_{n} 2√(N_{++,n} N_{−−,n}); ε_{n}: ułamek par tła wygrywających w klasie n. **Zielone leżą niżej niż czarne,** bo para prawdziwa, która przegrała ranking, już nie wraca.', { fontSize: 13 })], { x: 7.2, y: 4.82, w: 5.53, h: 1.2, fontSize: 13 });
  result(s, `zwykły wzór odejmuje o ${cbTooBig} % tła za dużo, a tła jest ${(nv(d.Ce)/nv(d.Se)).toFixed(1).replace('.', ',')} razy więcej niż sygnału.`);
}

{
  const pr = T.proj['0'];
  const s = cs(P5('Inne rozkłady'), 'W innych rozkładach k i wagi klas trzeba liczyć w tej samej zmiennej, bo k zmienia się nawet o 30 %',
    `Dotychczas k i ε_n były funkcjami masy pary. Tutaj sprawdzono, jak wyglądają w funkcji czterech innych zmiennych, dla par o masie 0,14–0,30 GeV/c², w danych PT3: masa brakująca pary e⁺e⁻, kąt pary w układzie środka masy, kąt e⁺ w układzie pary i pęd poprzeczny pary. Górny rząd to k z mieszania zdarzeń. Dolny rząd to ułamek par z isBest == 1 w zdarzeniach z trzema i czterema leptonami: pary w tych klasach to głównie tło, więc jest to bezpośredni, niezależny od symulacji obraz wag ε_n. Linie kropkowane to wartości średnie. k zmienia się bardzo: w kącie pary od ${pr.thcms.k[0]} do ${pr.thcms.k[1]}, w kącie e⁺ w układzie pary od ${pr.cosh.k[0]} do ${pr.cosh.k[1]}, w masie brakującej od ${pr.mm.k[0]} do ${pr.mm.k[1]}, a w p_T tylko od ${pr.pt.k[0]} do ${pr.pt.k[1]}. Skrajne wartości są na brzegach akceptancji. ε_n jest prawie stałe w masie brakującej i p_T, ale zmienia się z kątem pary, od ${pr.thcms.e3[0]} do ${pr.thcms.e3[1]} dla n = 3, bo ranking faworyzuje ślady w RPC. Dlatego wartości k(M) i ε_n(M) nie wolno przenosić na inne rozkłady: oba trzeba liczyć w zmiennej, która jest rysowana, z tymi samymi cięciami. Makro: plots/k_eps_projections.C.`);
  fig(s, 'figs/projections.png', MX, 1.5, 9.1, 4.2, 'Górny rząd: k z mieszania zdarzeń w funkcji masy brakującej, kąta pary, kąta e plus w układzie pary i pędu poprzecznego. Dolny rząd: ułamek par z isBest równym 1 dla trzech i czterech leptonów w tych samych zmiennych.', { left: true });
  cap(s, `Dane PT3, masa pary 0,14–0,30 ${M}. Dolny rząd: ułamek isBest == 1 w zdarzeniach z n = 3 i 4 leptonami. Linie kropkowane: średnie.`, MX, 5.55, 9.1, 0.45);
  bl(s, [
    '**Góra:** poprawka k w czterech innych zmiennych.',
    '**Dół:** jaka część par wygrywa ranking.',
    '**k** silnie zależy od kątów i masy brakującej, a słabo od p_{T}.',
    '**ε_{n}** jest prawie stałe poza zależnością od kąta pary.',
  ], { x: 9.95, y: 1.55, w: 2.78, h: 4.4, fontSize: 14, paraSpaceAfter: 6 });
  result(s, 'k i ε_{n} liczymy w tej zmiennej, którą rysujemy; wartości z widma masy nie wolno przenosić.');
}

// ================================================================ wnioski
section('Wnioski');
{
  const s = cs('Główny wniosek', 'Główny wniosek: isBest == 1 to para o najlepiej zgodnych czasach, a nie para najlepsza fizycznie',
    `To najważniejsza konkluzja studium. Znacznik isBest == 1 dostaje para o najlepiej zgodnych czasach przelotu, a kryterium nie zawiera informacji o tym, czy dwa leptony pochodzą z jednego rozpadu. Jeśli w zdarzeniu ranking wybierze parę tła, para prawdziwa wypada z próbki, a poprawne odjęcie tła jej nie przywróci. Korekta tła usuwa tylko tło, a ubytku sygnału nie naprawia. Dlatego zielone punkty na slajdzie 36 leżą niżej niż czarne. Strata dotyczy tylko par prawdziwych w zdarzeniach z co najmniej trzema leptonami, a tam jest około ${(100 - parseFloat(T.cls[0].sig.replace(',', '.'))).toFixed(0)} % sygnału. Pary w zdarzeniach z dwoma leptonami (${parseFloat(T.cls[0].sig.replace(',', '.')).toFixed(0)} % sygnału) nie mają rywala i zostają wszystkie. Tabela pokazuje dwa oszacowania tego, ile sygnału zostaje: pomiar to sygnał przy isBest == 1 z wagami i poprawką k podzielony przez sygnał przy isBest ≥ 0, a kolumna oczekiwana zakłada, że para prawdziwa zostaje tak często jak para tła. Niepewność wag to około 2 % sygnału, więc różnice między kolumnami nie są znaczące. Co z tym zrobić: do widma inkluzywnego brać wszystkie kombinacje, bo para prawdziwa zawsze zostaje w próbce. Gdy isBest == 1 jest konieczne, sygnał w klasie n można poprawić tylko średnio, dzieląc go przez ε_n. Gdy potrzebna jest właściwa para, trzeba ją wybierać kryterium fizycznym, na przykład masą e⁺e⁻γ albo masą brakującą.`);
  const lw = 7.15, g2 = 0.12; let y_ = 1.55;
  const hs = [1.08, 1.72, 1.32];
  card(s, 'Co naprawdę wybiera isBest == 1', 'Parę o najlepiej zgodnych czasach przelotu, a nie parę najlepszą fizycznie. Kryterium nie wie, czy leptony pochodzą z jednego rozpadu.', MX, y_, lw, hs[0], { size: 14.5, headSize: 16.5, fill: C.accent2 }); y_ += hs[0] + g2;
  card(s, 'Odjęcie tła nie przywraca utraconej pary prawdziwej', 'Jeśli w zdarzeniu ranking wybierze parę tła, para prawdziwa wypada z próbki, a poprawne odjęcie tła jej nie przywróci. Korekta tła usuwa tylko tło, a ubytku sygnału nie naprawia. Dlatego zielone punkty na slajdzie 36 leżą niżej niż czarne.', MX, y_, lw, hs[1], { size: 14.5, headSize: 16.5, fill: C.accent2 }); y_ += hs[1] + g2;
  card(s, 'Kogo dotyczy strata', `Tylko par prawdziwych w zdarzeniach z n ≥ 3, a tam jest około ${(100 - parseFloat(T.cls[0].sig.replace(',', '.'))).toFixed(0)} % sygnału. Pary w zdarzeniach z n = 2 (${parseFloat(T.cls[0].sig.replace(',', '.')).toFixed(0)} % sygnału) nie mają rywala i zostają wszystkie.`, MX, y_, lw, hs[2], { size: 14.5, headSize: 16.5 });
  const rx = MX + lw + 0.25, rw = CW - lw - 0.25;
  const rows = [['Przedział\n[GeV/c²]', { t: 'Zostaje\nsygnału', ...num }, { t: 'Oczekiwane', ...num }]];
  T.keep.forEach(r => rows.push([r.w, { t: r.meas, ...num, bold: true }, { t: r.exp, ...num }]));
  tbl(s, rows, { x: rx, y: 1.55, w: rw, colW: [1.75, 1.45, rw - 3.2], rowH: [0.56, 0.3, 0.3, 0.3, 0.3, 0.3], fontSize: 13, ruleRows: [4], boldRows: [5] });
  card(s, 'Co z tym zrobić', [
    '**Widmo inkluzywne:** isBest ≥ 0, para prawdziwa zawsze zostaje.',
    '**Gdy isBest == 1 jest konieczne:** sygnał w klasie n dzielić przez ε_{n}, tylko średnio.',
    '**Gdy potrzebna jest właściwa para:** kryterium fizyczne, np. masa e⁺e⁻γ.',
  ], rx, 3.98, rw, 2.0, { size: 13, headSize: 16.5, gap: 3, fill: C.accent1 });
  result(s, 'isBest ≥ 0 nie traci par prawdziwych; isBest == 1 traci 2–3 % sygnału i korekta tła tego nie naprawia.');
}
{
  const s = cs('Wnioski', 'Wnioski: wszystkie kombinacje do widma inkluzywnego; isBest == 1 traci sygnał i wymaga poprawki tła',
    'Do widma inkluzywnego należy brać wszystkie kombinacje: wzór jest wtedy poprawny do około 6 %, a resztę opisuje poprawka k. Gdy potrzebna jest dokładnie jedna para na zdarzenie, na przykład do rekonstrukcji π⁰, η i ω z fotonem w kalorymetrze, isBest == 1 pozwala wziąć jedną parę i zachowuje 99,6 % par prawdziwych, choć w zdarzeniach z kilkoma parami para prawdziwa nie zostaje częściej niż para tła. Tło trzeba wtedy liczyć osobno w klasach krotności, z wagami odtwarzającymi rywalizację. W podpróbce z trygerem PT3 taka poprawka zmienia stosunek wzoru do tła prawdziwego z 1,207 na 0,926, gdy wagi pochodzą z historii śladów, i na 0,918, gdy pochodzą z emulacji. Przy najmniejszych masach potrzebne jest cięcie na kąt otwarcia. Symulacja pokazuje mechanizm isBest, ale nie odtwarza tła kombinatorycznego w danych. Poprawkę k wyznaczono wyłącznie z mieszania zdarzeń w danych; wyniki w częściach 4 i 5. Kryterium wyboru isBest w kodzie PAT3 porównuje tylko czasy przelotu, więc nie może wskazać pary prawdziwej.');
  const w = 3.9, g = (CW - 3*w)/2, y = 1.55, h = 2.5;
  card(s, 'Widmo inkluzywne e⁺e⁻', `Wszystkie kombinacje (isBest ≥ 0). Wzór mnożymy przez poprawkę k z mieszania zdarzeń w danych: k = ${T.k_int.k} ± ${T.k_int.err} w 0,14–0,80 ${M}.`, MX, y, w, h, { fill: C.accent1 });
  card(s, 'Jedna para na zdarzenie', 'isBest == 1 to para o najlepiej zgodnych czasach, a nie najlepsza fizycznie. Tło trzeba liczyć z wagami ε_{n}, a para prawdziwa, która przegrała ranking, jest stracona: 2–3 % sygnału w danych w 0,14–0,80 GeV/c².', MX + w + g, y, w, h, { fill: C.accent2 });
  card(s, 'Najmniejsze masy', `Poniżej 0,02 ${M} przydaje się cięcie OA > 4°: stosunek wzoru do tła prawdziwego spada z 1,23 do 1,03 ± 0,04.`, MX + 2*(w + g), y, w, h);
  tx(s, '**Czego to studium nie rozstrzyga**', { x: MX, y: 4.28, w: CW, h: 0.32, fontSize: 17 });
  bl(s, [
    `Rola symulacji: pokazuje wpływ isBest, ale tła w danych nie odtwarza, więc k pochodzi z samych danych. Sam test metody w symulacji mówi, że mieszanie może nie obejmować całego k w 0,14–0,20 ${M}.`,
    'Sygnał w symulacji: isBest == 1 czy isBest ≥ 0 z pytaniem o kanał produkcji i rozpadu; do ustalenia.',
  ], { x: MX, y: 4.66, w: CW, h: 1.3 });
  result(s, 'isBest ≥ 0 nie traci par prawdziwych; isBest == 1 wymaga poprawki tła i traci 2–3 % sygnału.');
}

{
  const s = cs('Rekomendacje', 'Rekomendacje: jak analizować i jakie poprawki stosować',
    'Slajd zbiera zalecenia. isBest == 1 to para o najlepiej zgodnych czasach, a nie para najlepsza fizycznie, a para prawdziwa, która przegrała ranking, jest stracona. Do widma inkluzywnego należy brać wszystkie kombinacje i tło ze wzoru pomnożone przez k z mieszania zdarzeń w danych, przedział po przedziale. Gdy potrzebna jest jedna para na zdarzenie, isBest == 1 jest dopuszczalne, ale tło trzeba liczyć z wagami ε_n; zwykły wzór na parach isBest == 1 zaniża sygnał. Symulacja dostarcza tylko sygnału: pary prawdziwe wybrane po kanale produkcji i rozpadu, w tym samym wariancie isBest co dane. Tła kombinatorycznego z symulacji nie bierzemy. Gdy rysujemy inną zmienną niż masa albo stosujemy cięcie na kąt otwarcia, k i ε_n trzeba policzyć w tej zmiennej i z tym samym cięciem. Poprawka na tryger jest osobnym, kolejnym krokiem: mnoży gotowe widmo sygnału. Na tło, k i ε_n tryger nie wpływa: PT2 i PT3 dają te same wartości w granicach podanych na slajdzie.');
  tbl(s, [
    ['Zadanie', 'Pary', 'Tło', 'Skąd poprawka'],
    ['Widmo inkluzywne e⁺e⁻', 'isBest ≥ 0: para prawdziwa zawsze zostaje', 'k(M) · 2√(N_{++} N_{−−}), w każdym przedziale masy', 'k: mieszanie zdarzeń w danych'],
    ['Jedna para na zdarzenie, np. π⁰, η, ω → e⁺e⁻γ', 'kryterium fizyczne (masa e⁺e⁻γ, masa brakująca); isBest == 1 tylko z poprawkami', 'k · Σ_{n} ε_{n} 2√(N_{++,n} N_{−−,n}); nigdy zwykły wzór', 'k i ε_{n} z danych; stratę sygnału w klasie n poprawiać przez 1/ε_{n}'],
    ['Symulacja', 'pary prawdziwe po kanale produkcji i rozpadu, ten sam wariant isBest co w danych', 'tła z symulacji nie bierzemy', 'do porównań potrzebne przekroje czynne'],
    ['Inna zmienna niż masa', 'jak wyżej', 'k(X) i ε_{n}(X) policzone w tej zmiennej', 'te same pule i ranking, wartości z masy nie przenosić'],
    ['Cięcie na kąt otwarcia', 'jak wyżej', 'k liczone z tym samym cięciem', 'np. OA > 9°: k = 0,82 przy 10–30 MeV/c^{2}'],
  ], { x: MX, y: 1.5, w: CW, colW: [2.55, 3.2, 3.4, 2.98], rowH: [0.34, 0.5, 0.72, 0.72, 0.5, 0.5], fontSize: 12.5 });
  s.addText([...rt('**Tryger to osobny, następny krok.** Poprawka na tryger mnoży gotowe widmo sygnału. Na tło, k i ε_{n} tryger nie wpływa: PT2 i PT3 różnią się w k najwyżej o 0,4 % (poniżej 0,56 GeV/c^{2}), a w ε_{n} o 1,3 % dla n = 3 i 3 % dla n = 4.', { fontSize: 16 })],
    { x: MX, y: 5.2, w: CW, h: 0.86, color: C.text1, fill: { color: C.accent2, transparency: 75 }, margin: [6, 10, 6, 10], valign: 'middle', align: 'left', isTextBox: true, objectName: 'Tryger' });
  result(s, 'wszystkie kombinacje i k z danych; isBest == 1 tylko z wagami ε_{n}, ze świadomością straty sygnału.');
}

section('Dodatki');
{
  const s = cs('Dodatek', 'Tabela zapasowa: liczby par w przedziałach masy',
    'Tabela zapasowa zawiera liczby dla wszystkich przedziałów masy. Dla każdego przedziału podano liczbę wszystkich par e⁺e⁻, sygnał prawdziwy, tło prawdziwe, tło ze wzoru, ich stosunek z niepewnością statystyczną oraz stosunek sygnału ze wzoru do sygnału prawdziwego. Bez cięcia na kąt otwarcia.');
  const head = (t) => [t, { t: 'N_{+−}', ...num }, { t: 'Sygnał', ...num }, { t: 'Tło prawdziwe', ...num }, { t: 'Wzór', ...num }, { t: 'Wzór / prawdziwe', ...num }, { t: 'S_{wzór} / S', ...num }];
  const mk = (v, title) => { const rows = [head(title)]; T.windows.forEach(r => rows.push([r.w, { t: r[v].N, ...num }, { t: r[v].S, ...num }, { t: r[v].B, ...num }, { t: r[v].CB, ...num }, { t: `${r[v].R} ± ${r[v].eR}`, ...num }, { t: r[v].Sr, ...num }])); return rows; };
  const cw = [2.73, 1.6, 1.6, 1.6, 1.3, 1.9, 1.4];
  tbl(s, mk('all', 'isBest ≥ 0'), { x: MX, y: 1.5, w: CW, colW: cw, rowH: 0.3, fontSize: 13, boldRows: [6], ruleRows: [5] });
  tbl(s, mk('best', 'isBest == 1'), { x: MX, y: 3.85, w: CW, colW: cw, rowH: 0.3, fontSize: 13, boldRows: [6], ruleRows: [5] });
  cap(s, `Przedziały masy w ${M}, bez cięcia na kąt otwarcia, wszystkie zdarzenia. S_{wzór}: wszystkie pary minus tło ze wzoru, S: sygnał prawdziwy.`, MX, 6.2, CW, 0.3);
}

{
  const s = cs('Dodatek', 'Skład par jednoimiennych i tła e⁺e⁻',
    'Tabela przeniesiona ze slajdu o parach jednoimiennych. Składy pokazano w przedziale 0,14–0,80 GeV/c² dla isBest ≥ 0. Skład tła e⁺e⁻ leży między składem e⁺e⁺ i e⁻e⁻. Elektrony wtórne, z rozproszenia Comptona i elektrony delta, powstają tylko jako e⁻.');
  const a = T.ls.pp, b = T.ls.mm, c = T.cbcomp;
  tbl(s, [
    ['Skład [%]', { t: 'e⁺e⁺', ...num }, { t: 'e⁻e⁻', ...num }, { t: 'tło e⁺e⁻', ...num }],
    ['Dalitz + Dalitz', { t: a.dd, ...num }, { t: b.dd, ...num }, { t: c.dd, ...num }],
    ['Dalitz + konwersja', { t: a.dc, ...num }, { t: b.dc, ...num }, { t: c.dc, ...num }],
    ['konwersja + konwersja', { t: a.cc, ...num }, { t: b.cc, ...num }, { t: c.cc, ...num }],
    ['wtórne i inne', { t: a.sec, ...num }, { t: b.sec, ...num }, { t: c.sec, ...num }],
    ['fałszywy ślad', { t: a.fake, ...num }, { t: b.fake, ...num }, { t: c.fake, ...num }],
  ], { x: MX, y: 1.55, w: 7.4, colW: [3.5, 1.3, 1.3, 1.3], rowH: 0.48, fontSize: 16 });
  cap(s, `Przedział 0,14–0,80 ${M}, isBest ≥ 0. Liczby par: e⁺e⁺ ${a.tot}, e⁻e⁻ ${b.tot}, tło e⁺e⁻ 98 722.`, MX, 4.55, 7.4, 0.4);
  tx(s, 'Elektrony wtórne (rozproszenie Comptona, elektrony delta) występują tylko jako e⁻, dlatego ich udział wynosi 13,1 % w parach e⁻e⁻ i 0,6 % w parach e⁺e⁺.', { x: 8.4, y: 1.6, w: 4.33, h: 2.6 });
}

{
  const s = cs('Dodatek', 'Słownik pojęć',
    'Słownik zbiera pojęcia używane w prezentacji.');
  const left = [
    ['Pojęcie', 'Znaczenie'],
    ['tło kombinatoryczne', 'pary leptonów niepochodzących z jednego rozpadu ani jednego procesu'],
    ['pary różnoimienne', 'pary e⁺e⁻'],
    ['pary jednoimienne', 'pary e⁺e⁺ i e⁻e⁻'],
    ['hipoteza', 'jedna kombinacja dwóch leptonów w zdarzeniu'],
    ['isBest', '1 dla hipotezy o najlepiej zgodnych czasach przelotu (nie najlepszej fizycznie), 0 dla pozostałych'],
    ['klasa n', 'liczba zidentyfikowanych leptonów w zdarzeniu'],
    ['rozpad Dalitza', 'rozpad mezonu na foton i parę e⁺e⁻'],
    ['χ² czasów', 'miara zgodności czasów przelotu z hipotezą pary'],
  ];
  const right = [
    ['Pojęcie', 'Znaczenie'],
    ['konwersja', 'zamiana fotonu na parę e⁺e⁻ w materiale detektora'],
    ['kąt otwarcia (OA)', 'kąt między pędami dwóch leptonów pary'],
    ['β_{new}', 'prędkość wyznaczona z czasu względnego'],
    ['SMASH, GEANT', 'generator zderzeń; symulacja przejścia cząstek przez detektor'],
    ['PAT, FAT', 'programy analizy: zapis hipotez par oraz obróbka ntupli'],
    ['poprawka k', 'stosunek par e⁺e⁻ do wzoru w zdarzeniach mieszanych; opisuje różnicę akceptancji'],
    ['k_{true}', 'tło prawdziwe ze śladów podzielone przez wzór; tylko w symulacji'],
    ['mieszanie zdarzeń', 'łączenie leptonów z różnych zdarzeń; brak par skorelowanych'],
  ];
  tbl(s, left, { x: MX, y: 1.5, w: 5.95, colW: [2.1, 3.85], rowH: 0.5, fontSize: 14 });
  tbl(s, right, { x: MX + 6.18, y: 1.5, w: 5.95, colW: [2.1, 3.85], rowH: 0.5, fontSize: 14 });
}

// pptxgenjs writes <a:pPr> again before every later run of a paragraph: invalid XML, and the last one (no bullet) wins.
async function fixParagraphs(file) {
  const JSZip = require('jszip');
  const zip = await JSZip.loadAsync(fs.readFileSync(file));
  for (const name of Object.keys(zip.files).filter(n => /^ppt\/slides\/[^/]+\.xml$/.test(n))) {
    const xml = await zip.file(name).async('string');
    zip.file(name, xml.replace(/<a:p>(.*?)<\/a:p>/gs, (m, body) => '<a:p>' + body.replace(/<a:pPr\b[^>]*?(?:\/>|>.*?<\/a:pPr>)/gs, (pp, off) => (off === 0 ? pp : '')) + '</a:p>'));
  }
  fs.writeFileSync(file, await zip.generateAsync({ type: 'nodebuffer', compression: 'DEFLATE' }));
}
flushFootnotes();
(async () => {
  await pres.writeFile({ fileName: path.join(__dirname, OUT) });
  await fixParagraphs(path.join(__dirname, OUT));
  if (process.env.SKILL) { const { applyTheme } = require(path.join(process.env.SKILL, 'scripts/apply_theme.js')); await applyTheme(path.join(__dirname, OUT), L.THEME); }
  console.log('written', OUT);
})().catch(e => { console.error(e); process.exit(1); });
