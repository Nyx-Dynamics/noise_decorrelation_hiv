#!/usr/bin/env python3
"""
Locate run artifacts by partial timestamp or run name.

Default search root:
  quantum/quantum/results_v3_6/runs

Usage examples:
  python find_run_results.py 20251116_22284
  python find_run_results.py with_valcour_waic
  python find_run_results.py 20251117_1034 --json

Options:
  --runs-dir PATH   Override base runs directory
  --json            Emit JSON summary instead of human-readable text
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Any, List


def _default_runs_dir() -> Path:
    # Repo-root relative default
    here = Path(__file__).resolve().parent
    return here / 'quantum' / 'quantum' / 'results_v3_6' / 'runs'


def _collect_artifacts(run_dir: Path) -> Dict[str, Any]:
    info: Dict[str, Any] = {}
    run_info_path = run_dir / 'run_info.json'
    if run_info_path.exists():
        try:
            info.update(json.loads(run_info_path.read_text()))
        except Exception:
            pass

    # List key files
    traces = sorted([str(p) for p in run_dir.glob('trace*.nc')])
    summaries = sorted([str(p) for p in run_dir.glob('summary*.csv')])
    group_fits = sorted([str(p) for p in run_dir.glob('group_likelihood_fit*.csv')])
    ppc = sorted([str(p) for p in run_dir.glob('posterior_predictive*.csv')])
    results_ratio = sorted([str(p) for p in run_dir.glob('results_v3_6*_*.csv')])

    info.setdefault('paths', {})
    info['paths'] = {
        'run_dir': str(run_dir),
        'traces': traces,
        'summaries': summaries,
        'group_likelihood_fit': group_fits,
        'posterior_predictive': ppc,
        'other_results': results_ratio,
    }
    return info


def find_runs(query: str, runs_dir: Path) -> List[Dict[str, Any]]:
    results: List[Dict[str, Any]] = []
    if not runs_dir.exists():
        return results

    for sub in sorted(runs_dir.iterdir()):
        if not sub.is_dir():
            continue
        match = False
        run_info_path = sub / 'run_info.json'
        if query in sub.name:
            match = True
        ts = rn = None
        if run_info_path.exists():
            try:
                meta = json.loads(run_info_path.read_text())
                ts = str(meta.get('timestamp', ''))
                rn = str(meta.get('run_name', ''))
                if query in ts or query in rn:
                    match = True
            except Exception:
                pass
        if match:
            rec = _collect_artifacts(sub)
            if 'timestamp' not in rec and ts:
                rec['timestamp'] = ts
            if 'run_name' not in rec and rn:
                rec['run_name'] = rn
            results.append(rec)
    return results


def main(argv: List[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description='Locate run artifacts by partial timestamp or run name.')
    ap.add_argument('query', help='Substring to match (e.g., 20251116_22284 or with_valcour_waic)')
    ap.add_argument('--runs-dir', type=str, default=None, help='Override base runs directory')
    ap.add_argument('--json', action='store_true', help='Emit JSON output')
    args = ap.parse_args(argv)

    runs_dir = Path(args.runs_dir) if args.runs_dir else _default_runs_dir()
    matches = find_runs(args.query, runs_dir)

    if args.json:
        print(json.dumps({'query': args.query, 'runs_dir': str(runs_dir), 'matches': matches}, indent=2))
    else:
        if not matches:
            print(f"No runs matched '{args.query}' under {runs_dir}")
            return 1
        print(f"Found {len(matches)} match(es) for '{args.query}' under {runs_dir}:")
        for i, m in enumerate(matches, 1):
            rn = m.get('run_name') or Path(m['paths']['run_dir']).name
            ts = m.get('timestamp', '?')
            print(f"\n[{i}] {rn} @ {ts}")
            print(f"    run_dir: {m['paths']['run_dir']}")
            if m['paths']['traces']:
                print(f"    trace(s):")
                for p in m['paths']['traces']:
                    print(f"      - {p}")
            if m['paths']['summaries']:
                print(f"    summary CSV(s):")
                for p in m['paths']['summaries']:
                    print(f"      - {p}")
            if m['paths']['group_likelihood_fit']:
                print(f"    group fit CSV(s):")
                for p in m['paths']['group_likelihood_fit']:
                    print(f"      - {p}")
            if m['paths']['posterior_predictive']:
                print(f"    posterior predictive CSV(s):")
                for p in m['paths']['posterior_predictive']:
                    print(f"      - {p}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
