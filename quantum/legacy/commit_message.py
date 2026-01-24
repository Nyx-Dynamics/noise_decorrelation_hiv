#!/usr/bin/env python3
"""
Generate a Conventional Commits-style message from the current git working tree.

Usage:
  PYTHONPATH=. python Extra/commit_message.py [--type TYPE] [--scope SCOPE] [--summary TEXT] [--out PATH]

Notes:
- If --type/--scope/--summary are not provided, the script will infer a reasonable
  subject line from detected changes (Added/Modified/Deleted/Renamed) and paths.
- Writes the commit message to --out (default: results/commit_message.txt) and
  also prints it to stdout.

Examples:
  python Extra/commit_message.py
  python Extra/commit_message.py --type feat --scope bayes --summary "Add analyzer + generative viz utilities"
  python Extra/commit_message.py --out /tmp/commit_msg.txt
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from collections import defaultdict
from datetime import datetime

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

PORCELAIN_MAP = {
    'A': 'added',
    'M': 'modified',
    'D': 'deleted',
    'R': 'renamed',
    'C': 'copied',
    'U': 'updated',  # conflict/merge
    '?': 'untracked',
}

IMPORTANT_DIRS = (
    'quantum/',
    'Extra/',
    'datafiles/',
    'results/',
    'figures/',
    'project docs/',
)

IMPORTANT_FILES = (
    'Makefile',
    'requirements.txt',
    'requirements.lock.txt',
    'README.md',
    'project docs/venv_info.txt',
    'quantum/bayesian_optimization.py',
    'quantum/final_calibrated_model.py',
    'quantum/coupling_functions.py',
    'Extra/bayes_analyzer.py',
    'Extra/generative_viz.py',
)

SCOPE_HINTS = [
    ('quantum/bayesian_optimization.py', 'bayes'),
    ('Extra/bayes_analyzer.py', 'analysis'),
    ('Extra/generative_viz.py', 'viz'),
    ('quantum/', 'quantum'),
    ('Extra/', 'extra'),
    ('Makefile', 'build'),
    ('requirements', 'deps'),
    ('project docs', 'docs'),
]

def run(cmd: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


def parse_porcelain_z(output: str) -> list[dict]:
    # output is NUL-separated entries: XY NUL path [NUL path2] NUL
    parts = output.split('\0')
    entries = []
    i = 0
    while i < len(parts):
        rec = parts[i]
        i += 1
        if not rec:
            continue
        # rec starts with two status chars and a space before path
        if len(rec) < 3:
            continue
        x = rec[0]
        y = rec[1]
        rest = rec[3:] if len(rec) > 3 else ''
        path1 = rest
        path2 = None
        # Renames/Copies provide next field as new path
        if x in ('R', 'C') or y in ('R', 'C'):
            if i < len(parts):
                path2 = parts[i]
                i += 1
        entries.append({'index': x, 'worktree': y, 'path': path1, 'path2': path2})
    return entries


def collect_changes() -> tuple[dict, dict]:
    # Get branch
    branch = run(['git', 'rev-parse', '--abbrev-ref', 'HEAD']).stdout.strip() or '(detached)'

    # Get porcelain status (staged vs unstaged from first and second column)
    cp = run(['git', 'status', '--porcelain=v1', '-z'])
    if cp.returncode != 0:
        raise RuntimeError(f"git status failed: {cp.stderr}")
    entries = parse_porcelain_z(cp.stdout)

    staged = defaultdict(list)
    unstaged = defaultdict(list)

    for e in entries:
        # Prefer index status if present, else worktree
        if e['index'] and e['index'] != ' ':
            kind = PORCELAIN_MAP.get(e['index'], 'updated')
            key = f"{kind}"
            staged[key].append((e['path'], e['path2']))
        if e['worktree'] and e['worktree'] != ' ':
            kind = PORCELAIN_MAP.get(e['worktree'], 'updated')
            key = f"{kind}"
            unstaged[key].append((e['path'], e['path2']))

    return {'branch': branch, 'staged': staged}, {'unstaged': unstaged}


def infer_type(staged: dict, unstaged: dict) -> str:
    def has(kind: str) -> bool:
        return kind in staged and staged[kind] or kind in unstaged and unstaged[kind]
    # Priority: fix if only modifications to code mention bug fixes (heuristic unavailable), fallback to feat if any added python/code, else docs/chore
    # Heuristic: if only docs/Makefile/requirements changed -> docs or chore
    all_paths = []
    for coll in (staged, unstaged):
        for lst in coll.values():
            all_paths.extend([p for p, _ in lst if p])
    if any(p.endswith('.py') for p in all_paths) and any(p for p in all_paths if 'Extra/bayes_analyzer.py' in p or 'Extra/generative_viz.py' in p or 'quantum/bayesian_optimization.py' in p):
        return 'feat'
    if has('added'):
        return 'feat'
    # Docs-only?
    if all(p.endswith('.md') or p.endswith('.txt') or p in ('Makefile', 'requirements.txt', 'requirements.lock.txt') for p in all_paths):
        if any(p.endswith('.md') or p.endswith('.txt') for p in all_paths):
            return 'docs'
        return 'chore'
    # Default
    return 'chore'


def infer_scope(paths: list[str]) -> str:
    scopes = []
    for hint, scope in SCOPE_HINTS:
        if any(hint in p for p in paths):
            scopes.append(scope)
    # deduplicate and join
    scopes = sorted(set(scopes))
    return ','.join(scopes) if scopes else ''


def make_subject(commit_type: str, scope: str, summary: str) -> str:
    if scope:
        return f"{commit_type}({scope}): {summary}"
    return f"{commit_type}: {summary}"


def default_summary(paths: list[str]) -> str:
    added = [p for p in paths if p]
    # Provide a short generic summary based on hints
    bits = []
    if any('quantum/bayesian_optimization.py' in p for p in added):
        bits.append('Bayesian inference module')
    if any('Extra/bayes_analyzer.py' in p for p in added):
        bits.append('Bayes analyzer')
    if any('Extra/generative_viz.py' in p for p in added):
        bits.append('generative viz')
    if any('project docs/venv_info.txt' in p for p in added):
        bits.append('docs updates')
    if any(p == 'Makefile' for p in added):
        bits.append('Makefile targets')
    if any(p == 'requirements.txt' for p in added):
        bits.append('deps pin')
    if not bits:
        return 'update repository files'
    return 'add ' + ', '.join(bits)


def categorize(paths: list[tuple[str, str|None]]) -> dict:
    cats = defaultdict(list)
    for p, p2 in paths:
        group = 'other'
        if p.startswith('quantum/'):
            group = 'quantum'
        elif p.startswith('Extra/'):
            group = 'extra'
        elif p.startswith('project docs/') or p.endswith('.md') or p.endswith('.txt'):
            group = 'docs'
        elif p in ('Makefile', 'requirements.txt', 'requirements.lock.txt'):
            group = 'build'
        elif p.startswith('datafiles/') or p.startswith('results/') or p.startswith('figures/'):
            group = 'artifacts'
        cats[group].append((p, p2))
    return cats


def format_file_list(title: str, items: list[tuple[str, str|None]]) -> str:
    lines = [title]
    for p, p2 in items:
        if p2:
            lines.append(f"- {p} -> {p2}")
        else:
            lines.append(f"- {p}")
    return '\n'.join(lines)


def build_message(commit_type: str | None, scope: str | None, summary: str | None, out_path: str) -> str:
    meta_staged, meta_unstaged = collect_changes()
    staged = meta_staged['staged']
    unstaged = meta_unstaged['unstaged']

    # Aggregate paths
    all_paths = []
    for coll in (staged, unstaged):
        for lst in coll.values():
            all_paths.extend([p for p, _ in lst if p])

    inferred_type = commit_type or infer_type(staged, unstaged)
    inferred_scope = scope or infer_scope(all_paths)
    inferred_summary = summary or default_summary(all_paths)

    subject = make_subject(inferred_type, inferred_scope, inferred_summary)

    # Body header
    ts = datetime.now().strftime('%Y-%m-%d %H:%M')
    branch = meta_staged['branch']
    body = [
        f"Branch: {branch}",
        f"Date:   {ts}",
        "",
        "Summary of changes:",
    ]

    # Stage categories
    if staged:
        body.append("")
        body.append("Staged changes:")
        for kind, entries in staged.items():
            if not entries:
                continue
            cats = categorize(entries)
            for g, items in cats.items():
                body.append(format_file_list(f"  {kind} ({g})", items))
    if unstaged:
        body.append("")
        body.append("Unstaged changes:")
        for kind, entries in unstaged.items():
            if not entries:
                continue
            cats = categorize(entries)
            for g, items in cats.items():
                body.append(format_file_list(f"  {kind} ({g})", items))

    # Footer hints
    body.extend([
        "",
        "Notes:",
        "- Edit the subject/body as needed before committing.",
        "- Consider staging only the intended files before committing (git add ...).",
        "- Commit with: git commit -F results/commit_message.txt",
    ])

    message = subject + "\n\n" + "\n".join(body) + "\n"

    # Ensure results dir
    os.makedirs(os.path.join(ROOT, 'results'), exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write(message)

    return message


def main():
    parser = argparse.ArgumentParser(description='Generate a Conventional Commits-style message from git status')
    parser.add_argument('--type', dest='ctype', help='Conventional commit type (feat, fix, docs, chore, refactor, perf, test, build)')
    parser.add_argument('--scope', dest='scope', help='Scope to include in subject, e.g., bayes,quantum')
    parser.add_argument('--summary', dest='summary', help='Short subject summary text')
    parser.add_argument('--out', dest='out', default=os.path.join(ROOT, 'results', 'commit_message.txt'))
    args = parser.parse_args()

    try:
        msg = build_message(args.ctype, args.scope, args.summary, args.out)
        print(msg)
        print(f"\nCommit message written to: {args.out}")
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
