// lib_pl.js -- copy of lib.js (SLIDES_STYLE.md): FOOT, title, abbreviation file, Polish typography and "Wynik:" changed.
const fs = require('fs');
const path = require('path');
const pptxgen = require('pptxgenjs');
const HERE = __dirname;
const DATA = JSON.parse(fs.readFileSync(path.join(HERE, 'data/deck_data.json')));
const THEME = { name: 'HADES white', headFontFace: 'Calibri', bodyFontFace: 'Calibri',
  colors: { dk1: '111827', lt1: 'FFFFFF', dk2: '4B5563', lt2: 'F1F3F5', accent1: '2A78D6', accent2: 'EB6834', accent3: '1BAF7A', accent4: 'EDA100', accent5: '6250D6', accent6: '008300', hlink: '2A78D6', folHlink: '6250D6' } };
const HEX = { ink: '111827', ink2: '4B5563', grid: 'E5E7EB', rule: '9CA3AF', tint: 'F1F3F5', s1: '2A78D6', s2: 'EB6834', s3: '1BAF7A', s4: 'EDA100', s5: '6250D6' };
const ABBR = require('./abbrev_pl.js');
const pres = new pptxgen();
pres.layout = 'LAYOUT_WIDE';            // 13.33 x 7.5 in
pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
pres.title = 'Tło kombinatoryczne e⁺e⁻ w symulacji SMASH: sygnał, pary jednoimienne i wybór isBest';
pres.subject = 'Sygnał i tło kombinatoryczne par e⁺e⁻ w p+p, symulacja SMASH standard (gen4)';
pres.author = 'Analiza FAT, HADES feb22'; pres.company = 'HADES';
const C = pres.SchemeColor;
const W = 13.33, H = 7.5, MX = 0.6, CW = W - 2*MX;
const FOOT = 'HADES feb22 p+p, T = 4,53 GeV  ·  symulacja SMASH  ·  tło kombinatoryczne e⁺e⁻';
pres.defineSlideMaster({ title: 'CONTENT', background: { color: C.background1 }, slideNumber: { x: W - MX - 0.7, y: 7.17, w: 0.7, h: 0.25, fontSize: 11, color: C.text2, align: 'right' },
  objects: [
    { placeholder: { options: { name: 'title', type: 'title', x: MX, y: 0.52, w: CW, h: 0.98, fontSize: 28, bold: true, color: C.text1, align: 'left', valign: 'top', margin: 0 }, text: '' } },
    { placeholder: { options: { name: 'crumb', type: 'body', x: MX, y: 0.2, w: CW, h: 0.3, fontSize: 13, color: C.accent1, align: 'left', valign: 'top', margin: 0 }, text: '' } },
    { text: { text: FOOT, options: { x: MX, y: 7.17, w: 8, h: 0.25, fontSize: 11, color: C.text2, margin: 0, align: 'left' } } },
  ] });
pres.defineSlideMaster({ title: 'TITLE', background: { color: C.background1 },
  objects: [
    { placeholder: { options: { name: 'title', type: 'title', x: 0.9, y: 2.0, w: W - 1.8, h: 1.9, fontSize: 40, bold: true, color: C.text1, align: 'left', valign: 'bottom', margin: 0 }, text: '' } },
    { placeholder: { options: { name: 'subtitle', type: 'body', x: 0.9, y: 4.1, w: W - 1.8, h: 1.2, fontSize: 22, color: C.text2, align: 'left', valign: 'top', margin: 0 }, text: '' } },
  ] });
pres.defineSlideMaster({ title: 'DIVIDER', background: { color: C.background1 }, slideNumber: { x: W - MX - 0.7, y: 7.17, w: 0.7, h: 0.25, fontSize: 11, color: C.text2, align: 'right' },
  objects: [
    { placeholder: { options: { name: 'title', type: 'title', x: MX, y: 1.0, w: CW, h: 1.0, fontSize: 36, bold: true, color: C.text1, align: 'left', valign: 'bottom', margin: 0 }, text: '' } },
    { text: { text: FOOT, options: { x: MX, y: 7.17, w: 8, h: 0.25, fontSize: 11, color: C.text2, margin: 0, align: 'left' } } },
  ] });

// ---- rich text: **bold**, _{sub}, ^{sup}
// typesetting: a number stays with its unit, a label with its number, the last two words of a text stay together
const NBSP = ' ';
function typo(s) {
  return s.replace(/(\d(?:\}|\*\*)?) (%|nb|mb|μb|ns|MeV|GeV|mm|cm|m\b|kHz|Hz|mln|tys\.|par\b|plików|dni|sektorów|zdarzeń|śladów|leptonów|razy)/g, '$1' + NBSP + '$2')
    .replace(/ == /g, NBSP + '==' + NBSP).replace(/ ([±=<>≤≥≈∓]) /g, NBSP + '$1' + NBSP).replace(/(\d) × (\d)/g, '$1' + NBSP + '×' + NBSP + '$2')
    .replace(/(?<![\p{L}])(dzień|dni|sektor|sektorów|Część|slajd|slajdy|pasmo|plik|plików|klasa|klasy|kod|OA) (\d|>|I\b|II|III|IV\b)/gu, '$1' + NBSP + '$2')
    .replace(/(\d) (\d{3})(?!\d)/g, '$1' + NBSP + '$2').replace(/(\d) (\d{3})(?!\d)/g, '$1' + NBSP + '$2')
    .replace(/(?<=^|[\s\u00a0(„*])([wzioauWZIOAU]) /gu, '$1' + NBSP)
    .replace(/(\d)–(\d)/g, '$1\u2060–\u2060$2')
    .replace(/ (\S{1,14})$/, NBSP + '$1');
}
function rt(s0, base = {}) {
  const s = typo(s0);
  const out = [];
  const re = /(\*\*.+?\*\*|_\{.+?\}|\^\{.+?\})/g;
  let i = 0, m;
  const push = (t, o) => { if (t) out.push({ text: t, options: { ...base, ...o } }); };
  while ((m = re.exec(s)) !== null) {
    push(s.slice(i, m.index), {});
    const tok = m[0];
    if (tok.startsWith('**')) out.push(...rt(tok.slice(2, -2), { ...base, bold: true }));
    else if (tok.startsWith('_{')) push(tok.slice(2, -1), { subscript: true });
    else push(tok.slice(2, -1), { superscript: true });
    i = m.index + tok.length;
  }
  push(s.slice(i), {});
  return out.length ? out : [{ text: '', options: { ...base } }];
}
function bullets(items, base = {}) {
  const out = [];
  items.forEach((it, k) => {
    const sub = typeof it === 'object';
    const runs = rt(sub ? it.t : it, base);
    runs[0].options = { ...runs[0].options, bullet: sub ? { indent: 14 } : { indent: 15 }, indentLevel: sub ? 1 : 0 };
    if (k < items.length - 1) runs[runs.length - 1].options = { ...runs[runs.length - 1].options, breakLine: true };
    out.push(...runs);
  });
  return out;
}
function paras(items, base = {}) {
  const out = [];
  items.forEach((it, k) => {
    const runs = rt(it, base);
    if (k < items.length - 1) runs[runs.length - 1].options = { ...runs[runs.length - 1].options, breakLine: true };
    out.push(...runs);
  });
  return out;
}
function tx(slide, content, o) {
  const runs = typeof content === 'string' ? rt(content) : content;
  slide.addText(runs, { fontSize: 16, color: C.text1, margin: 0, valign: 'top', align: 'left', isTextBox: true, objectName: o.name || 'Text', ...o });
}
function bl(slide, items, o) { tx(slide, bullets(items), { paraSpaceAfter: 7, ...o, name: o.name || 'Bullets' }); }

// ---- abbreviations: footnote on the slide where they first appear (abbrev.js), also appended to the notes
const pendingFoot = [], usedAbbr = new Set();
const plain = s => s.replace(/\*\*/g, '').replace(/_\{(.+?)\}/g, '_$1').replace(/\^\{(.+?)\}/g, '^$1');
function withAbbr(key, notes) { return ABBR[key] ? notes + '\n\nSkróty na tym slajdzie: ' + plain(ABBR[key]) : notes; }
function footLater(slide, key) { if (ABBR[key]) { usedAbbr.add(key); pendingFoot.push([slide, key]); } }
function flushFootnotes() {
  const unused = Object.keys(ABBR).filter(k => !usedAbbr.has(k));
  if (unused.length) throw new Error('abbrev_pl.js: no slide with the title: ' + unused.join(' | '));
  pendingFoot.forEach(([slide, key]) => {
    const lines = Math.ceil(plain(ABBR[key]).length/138), h = 0.2*lines + 0.03;
    tx(slide, rt(ABBR[key]), { x: MX, y: 7.15 - h, w: CW, h, fontSize: 12, color: C.text2, valign: 'bottom', name: 'Abbreviations' });
  });
}

// ---- slides
let curSection = null;
function section(t) { pres.addSection({ title: t }); curSection = t; }
function cs(crumb, title, notes) {
  const s = pres.addSlide({ masterName: 'CONTENT', sectionTitle: curSection });
  s.addText(rt(title), { placeholder: 'title' });
  s.addText(crumb, { placeholder: 'crumb' });
  s.addNotes(withAbbr(title, notes));
  footLater(s, title);
  return s;
}
// ---- images
function pngSize(f) { const b = fs.readFileSync(f); return [b.readUInt32BE(16), b.readUInt32BE(20)]; }
function fig(slide, file, x, y, w, h, alt, o = {}) {
  const p = path.join(HERE, file); const [pw, ph] = pngSize(p); const r = pw/ph;
  let iw = w, ih = w/r; if (ih > h) { ih = h; iw = h*r; }
  const ix = o.left ? x : x + (w - iw)/2; const iy = o.top === false ? y + (h - ih)/2 : y;
  slide.addImage({ path: p, x: ix, y: iy, w: iw, h: ih, altText: alt, objectName: 'Figure: ' + file });
  return { x: ix, y: iy, w: iw, h: ih };
}
function cap(slide, t, x, y, w, h = 0.5) { tx(slide, t, { x, y, w, h, fontSize: 12, color: C.text2, name: 'Caption' }); }
const EQTXT = {};
function eq(slide, key, x, y, scale, alt, o = {}) {
  const [w0, h0] = DATA.eq_sizes[key]; const w = w0*scale, h = h0*scale;
  slide.addImage({ path: path.join(HERE, 'eq', key + '.png'), x: o.center ? x - w/2 : x, y, w, h, altText: alt, objectName: 'Equation ' + key });
  return { w, h };
}
function result(slide, s, y = 6.08, h = 0.52) {
  slide.addText(rt('**Wynik:** ' + s), { x: MX, y, w: CW, h, fontSize: 17, color: C.text1, fill: { color: C.background2 }, margin: [0, 10, 0, 10], valign: 'middle', align: 'left', isTextBox: true, objectName: 'Result' });
}
function line(s, base) { const r = rt(s, base); r[r.length - 1].options = { ...r[r.length - 1].options, breakLine: true }; return r; }   // one paragraph, however many runs the markup gives
function card(slide, head, body, x, y, w, h, o = {}) {
  const runs = [...line(head, { bold: true, fontSize: o.headSize || 17 }), ...(body ? paras(Array.isArray(body) ? body : [body], { fontSize: o.size || 15 }) : [])];
  slide.addText(runs, { x, y, w, h, color: C.text1, fill: o.fill ? { color: o.fill, transparency: 80 } : { color: C.background2 }, margin: [8, 10, 8, 10], valign: 'top', align: 'left', paraSpaceAfter: o.gap || 4, isTextBox: true, objectName: 'Card' });
}
function stat(slide, big, small, x, y, w, h = 1.25, o = {}) {
  slide.addText([...line(big, { fontSize: o.bigSize || 26, bold: true, color: o.color || C.text1 }), ...rt(small, { fontSize: o.size || 14, color: C.text2 })], { x, y, w, h, color: C.text1, fill: o.fill ? { color: o.fill, transparency: 80 } : { color: C.background2 }, margin: [6, 10, 6, 10], valign: 'middle', align: 'left', isTextBox: true, objectName: 'Stat' });
}
// ---- tables: rows = arrays of strings; first row is the header
function tbl(slide, rows, o) {
  const fs0 = o.fontSize || 14, aligns = o.aligns || [];
  const body = rows.map((r, i) => r.map((c, j) => {
    const cell = (typeof c === 'object' && !Array.isArray(c)) ? c : { t: c };
    const isH = i === 0 && !o.noHeader;
    const emph = (o.boldRows || []).includes(i);
    const sep = (o.ruleRows || []).includes(i);
    const noneB = { type: 'none' };
    const b = isH ? { type: 'solid', pt: 1, color: HEX.ink2 } : (sep ? { type: 'solid', pt: 0.75, color: HEX.rule } : { type: 'solid', pt: 0.25, color: HEX.grid });
    return { text: rt(String(cell.t), { fontSize: fs0, bold: isH || emph || cell.bold || undefined, color: isH ? C.text2 : C.text1 }),
      options: { align: cell.align || aligns[j] || 'left', valign: 'middle', border: [noneB, noneB, b, noneB], fill: cell.fill ? { color: cell.fill } : undefined, colspan: cell.colspan, margin: o.margin || [2, 6, 2, 6] } };
  }));
  slide.addTable(body, { x: o.x, y: o.y, w: o.w, colW: o.colW, rowH: o.rowH, autoPage: false, objectName: o.name || 'Table' });
}
// ---- charts
function chartBase(o = {}) {
  return { catAxisLabelFontFace: '+mn-lt', valAxisLabelFontFace: '+mn-lt', legendFontFace: '+mn-lt', dataLabelFontFace: '+mn-lt', catAxisTitleFontFace: '+mn-lt', valAxisTitleFontFace: '+mn-lt',
    catAxisLabelFontSize: 12, valAxisLabelFontSize: 12, legendFontSize: 13, dataLabelFontSize: 12, catAxisLabelColor: HEX.ink2, valAxisLabelColor: HEX.ink2, legendColor: HEX.ink, dataLabelColor: HEX.ink,
    valGridLine: { color: HEX.grid, size: 0.75 }, catGridLine: { style: 'none' }, showTitle: false, catAxisLineShow: true, valAxisLineShow: false,
    showCatAxisTitle: !!o.catTitle, catAxisTitle: o.catTitle || '', catAxisTitleFontSize: 13, catAxisTitleColor: HEX.ink2,
    showValAxisTitle: !!o.valTitle, valAxisTitle: o.valTitle || '', valAxisTitleFontSize: 13, valAxisTitleColor: HEX.ink2 };
}
module.exports = { withAbbr, footLater, flushFootnotes, fs, path, pres, C, HEX, THEME, DATA, W, H, MX, CW, rt, bullets, paras, tx, bl, section, cs, fig, cap, eq, result, card, stat, tbl, chartBase };
