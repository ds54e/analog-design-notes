# Reproduce the passive-condition study

Download the [passive study 1.0 example](https://ds54e.github.io/analog-design-notes/downloads/passive-study-v1.0.zip),
or use `examples/passive-sensitivity` in the public repository. The example
contains compact observations, exact original circuit bodies, saved MOS raw
parameters, the six fixed passive scenarios, and graph/statistics/replay code.
The public APM model is obtained separately under its original terms.

The commands below were actually run in the new directory
`/tmp/adn-passive-clean-001`, using only the selected example files plus a fresh
public model clone and Python environment. Use new destinations for another run.
The checked host was Linux x86_64, Python 3.9.25 and preinstalled ngspice 47 at
`/usr/local/bin/ngspice`. This is a clean input/Python check on the same host,
not a clean-OS or all-platform qualification. SPICE is outside the site build.

## Obtain dependencies

```bash
git -c credential.helper= -c http.extraHeader= clone --depth 1 --branch v5.0.0 https://github.com/ds54e/analog-process-models.git /tmp/adn-passive-clean-001/apm-v5
python3 -m venv /tmp/adn-passive-clean-001/.venv
/tmp/adn-passive-clean-001/.venv/bin/python -m pip install --disable-pip-version-check --no-cache-dir --index-url https://pypi.org/simple -r /tmp/adn-passive-clean-001/example/requirements.txt
```

The clone's annotated tag resolves to pinned commit
`381517fda5107fabf98af7801d5a5103f38e230c`. Replay checks that identity and the
APM045/VTG model file hashes. Install ngspice 47 separately when needed using
its [official source and documentation](https://ngspice.sourceforge.io/).
No private repository or GitHub authentication is required by the example.

## Regenerate the equations, statistics and graphs

Working directory for these commands: `/tmp/adn-passive-clean-001`.

```bash
.venv/bin/python -B example/statistics.py --data /tmp/adn-passive-clean-001/example/data --output /tmp/adn-passive-clean-001/recomputed-statistics.json
MPLCONFIGDIR=/tmp/adn-passive-clean-001/mplconfig PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -B example/plot.py --data /tmp/adn-passive-clean-001/example/data --output /tmp/adn-passive-clean-001/regenerated-plots
```

All six selected SVGs and the complete statistics JSON matched the source
bytes. The statistics script also reconstructs all 16 centered derivatives
from the 68 gradient observations and checks the local KCL estimates using
published nominal OP data. It then computes finite-cohort RMS/MAE, functional
counts, the conditional linear crossing and the exact two-sign decomposition.
There is no fitted passive distribution or new population confidence interval.
Graph regeneration does not rerun the underlying 5,724 native observations.

## Replay a saved physical circuit under a fixed passive condition

The following two native commands were executed in the same clean directory:

```bash
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -B example/replay.py --apm /tmp/adn-passive-clean-001/apm-v5 --output /tmp/adn-passive-clean-001/replay-R75-001 --candidate R75 --index 1000 --scenario ratio+0.1
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -B example/replay.py --apm /tmp/adn-passive-clean-001/apm-v5 --output /tmp/adn-passive-clean-001/replay-R37-001 --candidate R37 --index 1000 --scenario common-0.1
```

| Replayed observation | Output (V) | Gain at 1 kHz (V/V) | Current (µA) |
|---|---:|---:|---:|
| R75 / 1000 / ratio +0.10 | 0.5382745126518996 | −5.540382317205924 | 72.80971040271251 |
| R37 / 1000 / common −0.10 | 0.5536417417371056 | −5.390136792626277 | 36.99759867829463 |

All three scalars matched the original serialized observations exactly.
Both circuits still fail the assessed specification. **Successful reproduction
does not turn a physical failure into a passing design.**

Replay reconstructs the changed Rd/Rs body from the fixed log scenario and
checks its exact original target-body hash. It uses the CSV's saved MOS UIDs,
geometry, DELVTO and ln(MULU0); no sampling, trim, rebias or recentering occurs.
Pinned APM checks actual MOS model/geometry/m/nf/raw readbacks around OP/AC.
The example checks DC terminal domains, KCL and frequency endpoints. It also
reads effective Rs/Rd back from the actual voltage drops and measured terminal
or supply currents, using the documented resistor-load topology. The original
private acquisition separately recorded native `@resistor[resistance]` values.

## Public evidence and limits

- [Scientific evidence and clean replay receipt](evidence.json)
- [Frozen transfer plan](transfer-plan.json)
- [Complete new passive-condition observations](data/transfer.csv)
- [Known MOS baselines and OP-derived slopes](data/known-baselines.csv)
- [Saved physical MOS parameters](data/realizations.csv)
- [Nominal finite boundaries](data/nominal-boundaries.csv)
- [Data file manifest](data/manifest.json)

The MOS cohort is the already exposed first-study cohort, and the two-point
passive conditions are explicitly synthetic. Their finite counts do not estimate
manufacturing yield. The example replays selected OP/AC observations, not all
sine, noise, temperature or layout behavior. Full native observations remain
in the private experiment archive; they are not a public replay dependency.
The original [first study and its release](https://ds54e.github.io/analog-design-notes/studies/gain-operating-point-accuracy/)
remain available with their unchanged definitions and original download.
