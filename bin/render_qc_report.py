#!/usr/bin/env python3
"""
render_qc_report.py  —  Generate a per-assembly QC summary PDF.

Called by the QC_REPORT Nextflow process. All tool inputs are optional;
missing/sentinel files are silently skipped so the report degrades
gracefully when tools are disabled.

Pages
-----
1. Summary stats table (all parsed metrics)
2. kplexity curve + double-sigmoid fit  (only if kplex CSV is present)
"""

import argparse, sys, re
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit


# ── CLI ───────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--sample',                required=True)
    p.add_argument('--type',                  required=True, dest='asm_type')
    p.add_argument('--assembly_fasta',        required=True)
    p.add_argument('--outfile',               required=True)
    p.add_argument('--quast_tsv',             default=None)
    p.add_argument('--merqury_qv',            default=None)
    p.add_argument('--merqury_completeness',  default=None)
    p.add_argument('--compleasm_summary',     default=None)
    p.add_argument('--lai_output',            default=None)
    p.add_argument('--kplex_csv',             default=None)
    p.add_argument('--flagger_bed',           default=None)
    return p.parse_args()


def _is_real(path):
    """Return True if path is a real, non-empty file (not the NO_FILE sentinel)."""
    if not path:
        return False
    p = Path(path)
    return p.exists() and p.stat().st_size > 0 and p.name != 'NO_FILE'


# ── FASTA N50 ─────────────────────────────────────────────────────────────────
def fasta_stats(fa_path):
    """Return dict: n_seqs, n50_mb, l50, total_mb."""
    lengths = []
    cur = 0
    with open(fa_path) as fh:
        for line in fh:
            if line.startswith('>'):
                if cur:
                    lengths.append(cur)
                cur = 0
            else:
                cur += len(line.strip())
    if cur:
        lengths.append(cur)
    lengths.sort(reverse=True)
    total = sum(lengths)
    running = 0
    for i, l in enumerate(lengths):
        running += l
        if running >= total / 2:
            return dict(n_seqs=len(lengths),
                        n50_mb=round(l / 1e6, 2),
                        l50=i + 1,
                        total_mb=round(total / 1e6, 2))
    return dict(n_seqs=len(lengths), n50_mb=0, l50=0, total_mb=round(total / 1e6, 2))


# ── QUAST ─────────────────────────────────────────────────────────────────────
def parse_quast(tsv_path):
    """Parse QUAST report.tsv → dict of key metrics."""
    metrics = {}
    want = {
        '# contigs': 'n_contigs',
        'Total length': 'total_length_bp',
        'N50': 'n50_bp',
        'L50': 'l50',
        'N75': 'n75_bp',
        'Genome fraction (%)': 'genome_fraction_pct',
        'Duplication ratio': 'dup_ratio',
        '# misassemblies': 'misassemblies',
        '# mismatches per 100 kbp': 'mismatches_per_100k',
        '# indels per 100 kbp': 'indels_per_100k',
    }
    with open(tsv_path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            key = parts[0].strip()
            val = parts[1].strip()
            for pattern, out_key in want.items():
                if key.startswith(pattern):
                    try:
                        metrics[out_key] = float(val.replace(',', ''))
                    except ValueError:
                        metrics[out_key] = val
                    break
    return metrics


# ── Merqury ───────────────────────────────────────────────────────────────────
def parse_merqury_qv(qv_path):
    """Return (qv, error_rate) from the last data line of a Merqury .qv file."""
    last = None
    with open(qv_path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith('#'):
                last = line
    if not last:
        return None, None
    parts = last.split()
    try:
        return float(parts[-2]), float(parts[-1])
    except (IndexError, ValueError):
        return None, None


def parse_merqury_completeness(stats_path):
    """Return completeness % from .completeness.stats (last data line, col 5)."""
    last = None
    with open(stats_path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith('#'):
                last = line
    if not last:
        return None
    parts = last.split()
    try:
        return float(parts[4])
    except (IndexError, ValueError):
        return None


# ── Compleasm ─────────────────────────────────────────────────────────────────
def parse_compleasm(summary_path):
    """Parse compleasm summary.txt → dict with S, D, F, M percentages."""
    result = {}
    pct_re = re.compile(r'([SDFMI]):\s*(\d+),\s*([\d.]+)%')
    with open(summary_path) as fh:
        for line in fh:
            m = pct_re.search(line)
            if m:
                result[m.group(1)] = float(m.group(3))
    return result


# ── LAI ───────────────────────────────────────────────────────────────────────
def parse_lai(lai_path):
    """Return whole-genome LAI score (float) or None."""
    with open(lai_path) as fh:
        for line in fh:
            if 'whole_genome' in line:
                parts = line.split()
                try:
                    return float(parts[6])
                except (IndexError, ValueError):
                    pass
    return None


# ── Flagger ───────────────────────────────────────────────────────────────────
def parse_flagger_bed(bed_path):
    """
    Parse a Flagger BED file.
    Returns dict: label → percentage of total assembly length.
    Expects column 4 to carry annotation labels (Hap, Collapse, Error, Duplication, etc.).
    """
    totals = {}
    grand = 0
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('track'):
                continue
            parts = line.split('\t')
            if len(parts) < 4:
                continue
            try:
                length = int(parts[2]) - int(parts[1])
            except ValueError:
                continue
            label = parts[3].strip()
            totals[label] = totals.get(label, 0) + length
            grand += length
    if grand == 0:
        return {}
    return {k: round(100 * v / grand, 2) for k, v in totals.items()}


# ── kplexity double-sigmoid fit ───────────────────────────────────────────────
def _double_sigmoid(x, L1, s1, k1, L2, s2, k2):
    return L1 / (1 + np.exp(-s1 * (x - k1))) + L2 / (1 + np.exp(-s2 * (x - k2)))


def fit_kplexity(csv_path, n_restarts=20):
    """
    Read kplexity CSV, fit double-sigmoid, return (k_vals, y_vals, popt, asymptote, r2).
    Returns (None, None, None, None, None) on failure.
    """
    import csv as _csv
    rows = []
    with open(csv_path) as fh:
        reader = _csv.DictReader(fh)
        fieldnames = reader.fieldnames or []
        # find k and fraction_unique columns flexibly
        k_col = next((c for c in fieldnames if c.lower() in ('k', 'kmer', 'kmer_size', 'k_size')), None)
        y_col = next((c for c in fieldnames if 'unique' in c.lower() or 'fraction' in c.lower()), None)
        if k_col is None or y_col is None:
            # fallback: use first two numeric columns
            k_col = fieldnames[0] if fieldnames else None
            y_col = fieldnames[1] if len(fieldnames) > 1 else None
        if k_col is None:
            return None, None, None, None, None
        for row in reader:
            try:
                rows.append((float(row[k_col]), float(row[y_col])))
            except (KeyError, ValueError):
                pass

    if len(rows) < 10:
        return None, None, None, None, None

    k_vals = np.array([r[0] for r in rows])
    y_vals = np.array([r[1] for r in rows])

    best_popt, best_r2 = None, -np.inf
    rng = np.random.default_rng(42)

    for _ in range(n_restarts):
        p0 = [
            rng.uniform(0.3, 0.7),   # L1
            rng.uniform(0.05, 0.5),  # s1
            rng.uniform(5, 30),      # k0_1
            rng.uniform(0.1, 0.5),   # L2
            rng.uniform(0.05, 0.3),  # s2
            rng.uniform(40, 120),    # k0_2
        ]
        bounds = ([0,0,0, 0,0,0], [1,5,60, 1,5,151])
        try:
            popt, _ = curve_fit(_double_sigmoid, k_vals, y_vals,
                                p0=p0, bounds=bounds, maxfev=10000)
            y_pred = _double_sigmoid(k_vals, *popt)
            ss_res = np.sum((y_vals - y_pred) ** 2)
            ss_tot = np.sum((y_vals - np.mean(y_vals)) ** 2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
            if r2 > best_r2:
                best_r2 = r2
                best_popt = popt
        except Exception:
            continue

    if best_popt is None:
        return k_vals, y_vals, None, None, None

    asymptote = best_popt[0] + best_popt[3]
    return k_vals, y_vals, best_popt, round(asymptote, 4), round(best_r2, 4)


# ── PDF pages ─────────────────────────────────────────────────────────────────
def draw_stats_page(fig, sample, asm_type, data):
    """Render a two-column key/value table of all QC metrics."""
    ax = fig.add_axes([0.05, 0.05, 0.90, 0.88])
    ax.axis('off')
    ax.set_title(f'QC Summary — {sample}  [{asm_type}]',
                 fontsize=14, fontweight='bold', pad=16)

    # Build display rows
    rows = []
    SEP = ('─────────────────────', '')

    # Assembly size
    if 'fasta' in data:
        fs = data['fasta']
        rows += [
            ('Assembly', ''),
            ('  Total length (Mb)', f"{fs['total_mb']}"),
            ('  Sequences', str(fs['n_seqs'])),
            ('  N50 (Mb)', f"{fs['n50_mb']}"),
            ('  L50', str(fs['l50'])),
            SEP,
        ]

    # QUAST
    if 'quast' in data:
        q = data['quast']
        rows.append(('QUAST', ''))
        for label, key in [
            ('  Genome fraction (%)', 'genome_fraction_pct'),
            ('  Misassemblies', 'misassemblies'),
            ('  Mismatches / 100 kbp', 'mismatches_per_100k'),
            ('  Indels / 100 kbp', 'indels_per_100k'),
        ]:
            if key in q:
                rows.append((label, str(q[key])))
        rows.append(SEP)

    # Merqury
    if 'merqury_qv' in data or 'merqury_completeness' in data:
        rows.append(('Merqury', ''))
        if 'merqury_qv' in data:
            qv, err = data['merqury_qv']
            rows.append(('  QV', f"{qv:.1f}"))
            if err is not None:
                rows.append(('  Error rate', f"{err:.2e}"))
        if 'merqury_completeness' in data:
            rows.append(('  K-mer completeness (%)', f"{data['merqury_completeness']:.2f}"))
        rows.append(SEP)

    # Compleasm
    if 'compleasm' in data:
        c = data['compleasm']
        rows.append(('BUSCO (compleasm)', ''))
        for tag in ['S', 'D', 'F', 'M']:
            if tag in c:
                label_map = {'S': 'Single', 'D': 'Duplicated', 'F': 'Fragmented', 'M': 'Missing'}
                rows.append((f"  {label_map[tag]} (%)", f"{c[tag]:.1f}"))
        rows.append(SEP)

    # LAI
    if 'lai' in data:
        rows += [
            ('LAI (LTR Assembly Index)', ''),
            ('  Score', f"{data['lai']:.2f}"),
            ('  Grade', 'Reference' if data['lai'] >= 10 else 'Acceptable' if data['lai'] >= 5 else 'Draft'),
            SEP,
        ]

    # kplexity
    if 'kplex_asymptote' in data:
        rows += [
            ('kplexity (k=5–151)', ''),
            ('  Unique-k asymptote', f"{data['kplex_asymptote']:.4f}"),
            ('  Fit R²', f"{data['kplex_r2']:.4f}"),
            SEP,
        ]

    # Flagger
    if 'flagger' in data:
        rows.append(('Flagger', ''))
        for label, val in sorted(data['flagger'].items()):
            rows.append((f"  {label} (%)", f"{val:.2f}"))
        rows.append(SEP)

    if not rows:
        ax.text(0.5, 0.5, 'No QC data available.', ha='center', va='center', fontsize=12)
        return

    # Remove trailing separator
    while rows and rows[-1] == SEP:
        rows.pop()

    col_labels = ['Metric', 'Value']
    tbl = ax.table(cellText=rows, colLabels=col_labels,
                   loc='upper center', cellLoc='left',
                   colWidths=[0.7, 0.3])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 1.6)

    # Style header
    for col in range(2):
        tbl[(0, col)].set_facecolor('#2c3e50')
        tbl[(0, col)].set_text_props(color='white', fontweight='bold')

    # Style rows
    for i in range(1, len(rows) + 1):
        metric, val = rows[i - 1]
        is_sep = metric.startswith('─')
        is_header = (val == '' and not is_sep)
        for col in range(2):
            cell = tbl[(i, col)]
            if is_sep:
                cell.set_facecolor('#ecf0f1')
                cell.set_text_props(color='#7f8c8d')
            elif is_header:
                cell.set_facecolor('#d5e8f3')
                cell.set_text_props(fontweight='bold')
            else:
                cell.set_facecolor('#ffffff' if i % 2 == 0 else '#f8f9fa')
            cell.set_edgecolor('#bdc3c7')


def draw_kplex_page(fig, sample, asm_type, k_vals, y_vals, popt, asymptote, r2):
    """Render kplexity curve with double-sigmoid fit."""
    ax = fig.add_subplot(111)

    ax.scatter(k_vals, y_vals, s=12, color='#3498db', alpha=0.6,
               zorder=3, label='Observed')

    if popt is not None:
        k_smooth = np.linspace(k_vals.min(), k_vals.max(), 500)
        y_smooth = _double_sigmoid(k_smooth, *popt)
        ax.plot(k_smooth, y_smooth, color='#e74c3c', lw=2,
                label=f'Double-sigmoid fit\nasymptote = {asymptote:.4f}  R² = {r2:.4f}')
        ax.axhline(asymptote, color='#e74c3c', lw=1, ls='--', alpha=0.5)
        ax.text(k_vals.max() * 0.98, asymptote + 0.005,
                f'asymptote = {asymptote:.4f}',
                ha='right', va='bottom', fontsize=9, color='#c0392b')

    ax.set_xlabel('k-mer size (k)', fontsize=11)
    ax.set_ylabel('Fraction unique k-mers', fontsize=11)
    ax.set_title(f'kplexity curve — {sample}  [{asm_type}]',
                 fontsize=13, fontweight='bold')
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=9, loc='lower right')
    ax.grid(True, alpha=0.3)
    fig.tight_layout()


# ── Main ─────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()

    data = {}

    # FASTA stats (always available)
    try:
        data['fasta'] = fasta_stats(args.assembly_fasta)
    except Exception as e:
        print(f"[WARN] FASTA stats failed: {e}", file=sys.stderr)

    # QUAST
    if _is_real(args.quast_tsv):
        try:
            data['quast'] = parse_quast(args.quast_tsv)
        except Exception as e:
            print(f"[WARN] QUAST parse failed: {e}", file=sys.stderr)

    # Merqury QV
    if _is_real(args.merqury_qv):
        try:
            data['merqury_qv'] = parse_merqury_qv(args.merqury_qv)
        except Exception as e:
            print(f"[WARN] Merqury QV parse failed: {e}", file=sys.stderr)

    # Merqury completeness
    if _is_real(args.merqury_completeness):
        try:
            data['merqury_completeness'] = parse_merqury_completeness(args.merqury_completeness)
        except Exception as e:
            print(f"[WARN] Merqury completeness parse failed: {e}", file=sys.stderr)

    # Compleasm
    if _is_real(args.compleasm_summary):
        try:
            data['compleasm'] = parse_compleasm(args.compleasm_summary)
        except Exception as e:
            print(f"[WARN] Compleasm parse failed: {e}", file=sys.stderr)

    # LAI
    if _is_real(args.lai_output):
        try:
            score = parse_lai(args.lai_output)
            if score is not None:
                data['lai'] = score
        except Exception as e:
            print(f"[WARN] LAI parse failed: {e}", file=sys.stderr)

    # kplexity
    kplex_result = (None, None, None, None, None)
    if _is_real(args.kplex_csv):
        try:
            kplex_result = fit_kplexity(args.kplex_csv)
            k_vals, y_vals, popt, asymptote, r2 = kplex_result
            if asymptote is not None:
                data['kplex_asymptote'] = asymptote
                data['kplex_r2'] = r2
        except Exception as e:
            print(f"[WARN] kplexity fit failed: {e}", file=sys.stderr)

    # Flagger
    if _is_real(args.flagger_bed):
        try:
            data['flagger'] = parse_flagger_bed(args.flagger_bed)
        except Exception as e:
            print(f"[WARN] Flagger parse failed: {e}", file=sys.stderr)

    # ── Render PDF ────────────────────────────────────────────────────────────
    k_vals, y_vals, popt, asymptote, r2 = kplex_result
    has_kplex = k_vals is not None and len(k_vals) > 0

    with PdfPages(args.outfile) as pdf:
        # Page 1: stats table
        fig = plt.figure(figsize=(10, 12))
        draw_stats_page(fig, args.sample, args.asm_type, data)
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # Page 2: kplexity curve (only if data present)
        if has_kplex:
            fig = plt.figure(figsize=(10, 6))
            draw_kplex_page(fig, args.sample, args.asm_type,
                            k_vals, y_vals, popt, asymptote, r2)
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

        d = pdf.infodict()
        d['Title'] = f'QC Report — {args.sample} [{args.asm_type}]'
        d['Author'] = 'asm-qc-nf'

    print(f"[render_qc_report] Written: {args.outfile}")


if __name__ == '__main__':
    main()
