import glob
import os

"""
Tiny smoke test for geometry comparison utilities.

Usage:
  PYTHONPATH=. python Extra/geometry_compare_smoke_test.py

It finds the latest *_summary.json and *_sse_coherence.csv in datafiles/,
then runs the comparison (Markdown + JSON) and the viz geometry plot.
"""


def latest(pattern: str) -> str:
    files = glob.glob(pattern)
    if not files:
        return None
    files.sort(key=lambda p: os.path.getmtime(p), reverse=True)
    return files[0]


def run():
    summary = latest('datafiles/*_summary.json')
    if not summary:
        raise SystemExit('No *_summary.json found under datafiles/. Run a simulation first.')
    csvp = latest('datafiles/*_sse_coherence.csv')

    from quantum.geometry_compare import compare_single, plot_geometry

    compare_single(summary, csv_hint=csvp)

    if csvp:
        plot_geometry(csvp)
    else:
        print('No *_sse_coherence.csv found; skipped geometry plot.')

    print('Geometry comparison smoke test completed.')


if __name__ == '__main__':
    run()
