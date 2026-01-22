import glob
import json
import os

"""
Tiny smoke test for the data interpretation mechanism.

Usage:
  PYTHONPATH=. python Extra/sse_interpret_smoke_test.py

It finds the latest *_summary.json in datafiles/, runs the interpreter,
then checks that Markdown and JSON interpretation outputs were written.
"""

def latest_summary(pattern: str = 'datafiles/*_summary.json'):
    files = glob.glob(pattern)
    if not files:
        return None
    files.sort(key=lambda p: os.path.getmtime(p), reverse=True)
    return files[0]


def run():
    path = latest_summary()
    if not path:
        raise SystemExit('No *_summary.json found under datafiles/. Run a simulation first.')

    from quantum.sse_interpret import interpret_run
    findings = interpret_run(path)
    base = os.path.basename(path).replace('_summary.json', '')
    md = os.path.join('results', f'{base}_interpretation.md')
    js = os.path.join('results', f'{base}_interpretation.json')
    assert os.path.exists(md), f"Markdown interpretation not found: {md}"
    assert os.path.exists(js), f"JSON interpretation not found: {js}"
    with open(js, 'r') as f:
        data = json.load(f)
    # Basic keys present
    for k in ['hiv_phase', 'dephasing_model', 'gamma_Tegmark', 'stability_flag']:
        assert k in data, f"Missing key in interpretation JSON: {k}"
    print('Interpretation smoke test passed for', path)


if __name__ == '__main__':
    run()
