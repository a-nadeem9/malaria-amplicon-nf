"""Observed first/last clonal turnover, with explicit per-locus sampling bounds."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

from .longitudinal import build_observed_dynamics


def calculate_clonal_dynamics(calls: pd.DataFrame) -> tuple[dict[str, pd.DataFrame], dict]:
    required = {'participant_id', 'collection_date', 'locus', 'allele'}
    missing = required - set(calls.columns)
    if missing:
        raise ValueError('Clonal input is missing: ' + ', '.join(sorted(missing)))
    data = calls.copy()
    valid = pd.Series(True, index=data.index)
    for col in required:
        valid &= data[col].notna() & data[col].astype(str).str.strip().ne('')
    dates = pd.to_datetime(data.collection_date, errors='coerce')
    valid &= dates.notna()
    excluded = data.loc[~valid].copy()
    data = data.loc[valid].copy()
    data['collection_date'] = dates[valid].dt.strftime('%Y-%m-%d')
    if 'max_abundance_pct' not in data:
        data['max_abundance_pct'] = float('nan')
    tables: dict[str, list[pd.DataFrame]] = {}
    for locus, group in data.groupby('locus', sort=True):
        # A visit at a different locus cannot establish this locus's absence.
        for name, table in build_observed_dynamics(group).items():
            table = table.copy()
            table['locus'] = locus
            tables.setdefault(name, []).append(table)
    output = {name: pd.concat(parts, ignore_index=True) for name, parts in tables.items()}
    output['excluded_calls'] = excluded
    turnover = output.get('turnover_metrics', pd.DataFrame())
    eligible = turnover[turnover.eligible_longitudinal] if not turnover.empty else turnover
    summary = {
        'schema_version': 1,
        'method': 'Observed first/last appearance per participant and locus',
        'acquisition_rule': 'First allele observation is later than the first genotyped visit at this locus.',
        'clearance_rule': 'Last allele observation is earlier than the last genotyped visit at this locus.',
        'denominator': 'max(1, distinct genotyped dates at this locus minus 1)',
        'interpretation': 'Changes in observed genotypes, not proof of biological acquisition or clearance. Missing genotypes and PCR-only visits are excluded from these rates. Reappearance of the same allele does not create another first/last event.',
        'single_visit_rule': 'Historical formula returns zero; eligible_longitudinal is false and these rows are excluded from rate summaries.',
        'participants': int(data.participant_id.nunique()),
        'loci': sorted(map(str, data.locus.unique())),
        'participant_locus_rows': len(turnover),
        'eligible_participant_locus_rows': len(eligible),
        'excluded_calls': len(excluded),
        'acquisitions': int(eligible.acquisitions.sum()) if len(eligible) else 0,
        'clearances': int(eligible.clearances.sum()) if len(eligible) else 0,
        'median_turnover_per_interval': float(eligible.turnover_per_interval.median()) if len(eligible) else None,
    }
    return output, summary


def export_clonal_dynamics(source: Path, output_dir: Path, mode: str) -> dict:
    calls = pd.read_csv(source, sep='\t')
    tables, summary = calculate_clonal_dynamics(calls)
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, table in tables.items():
        table.to_csv(output_dir / f'{name}.tsv', sep='\t', index=False)
    summary.update({'analysis_mode': mode, 'source': str(source),
                    'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
                    'files': {name: str(output_dir / f'{name}.tsv') for name in tables}})
    (output_dir / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    preview = tables.get('turnover_metrics', pd.DataFrame()).head(100)
    return {'summary': summary, 'rows': json.loads(preview.to_json(orient='records'))}
