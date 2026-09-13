# Reproduce the design-decisions supplement

This **learning-decisions-v1.0** supplement accompanies learning edition **v1.2**.
Download [the selected ZIP](https://ds54e.github.io/analog-design-notes/downloads/learning-decisions-v1.0.zip), which
contains one `example/` directory, or inspect the
[public files](https://github.com/ds54e/analog-design-notes/tree/main/examples/learning-decisions-v1).
All seven earlier ZIPs and their publication tags remain unchanged.

The supply/loss reanalysis supports chapters 6 and 7. It uses **known observations**
from the original finite-feedback study. It draws no new physical samples and
does not change that study's confirmation, definitions, circuits or release.

The worked result has two parts. A four-node small-signal KCL calculation
separates finite-reference supply injection from its output-voltage path.
An empirical error calculation then explains why F75L can have lower SD,
RMS and MAE at 0.97 V while fewer of its known circuits meet the fixed center
limit than F75S. Neither result is a new statistical confirmation.

## Exact inputs

The required input is the preserved
[feedback-study-v1.0.zip](https://ds54e.github.io/analog-design-notes/downloads/feedback-study-v1.0.zip), SHA256
`348d3f2c8cb3ef41b9c968f4bbb68c011c40cfb64ddda4138d2ff115c5e95bc6`.
It belongs to tag `finite-feedback-v1.0`, public commit
`721a59b9c33146b6ff1d8662e6b9a5076aafe448`. The analysis reads fifteen explicit
members and checks their hashes against [the input record](data/source.json)
and the original data manifest. It does not extract or execute the old package.

Those selected members contain four exact circuit definitions, nominal OP and
metrics, local sensitivity and condition tables, the original contract and
measurement/confirmation definitions, and 4,168 known DC/AC/current condition
records. The original 521 coordinates, candidate inventory differences,
physical realization IDs, and failed performance checks remain intact.

The original models were APM v5.0.0/APM045 VTG with ngspice 47. This reanalysis
needs only saved observations; no simulator or model installation is needed.
Use the [original study reproduction](https://ds54e.github.io/analog-design-notes/studies/finite-feedback/reproduce.html) or
the separate [initial-OTA example](https://ds54e.github.io/analog-design-notes/studies/learning-day-v1/reproduce.html) for actual
native replay, which is a different operation from these calculations.

## Recalculate the supply/loss explanation and its two graphs

Run from a new directory with the selected example files. Use Python 3.9 or
later and the pinned requirements supplied with them. Retrieve the original
ZIP at the URL above; no GitHub account is required.

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r example/requirements.txt
curl --fail --location https://ds54e.github.io/analog-design-notes/downloads/feedback-study-v1.0.zip --output feedback-study-v1.0.zip
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B example/analyze.py --feedback-zip feedback-study-v1.0.zip --output "$PWD/analysis-001"
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B example/plot.py --feedback-zip feedback-study-v1.0.zip --output "$PWD/figures-001"
```

The scripts refuse existing output destinations. `analyze.py` writes
`supply-paths.csv`, `loss-summary.csv`, `paired-supply.csv`, `known-errors.csv`
and `summary.json`. `plot.py` produces SVG/PNG versions of
`supply-reference-paths` and `fixed-target-loss`, directly from the same
verified original ZIP. The latter figure is an empirical cumulative fraction;
there is no fitted normal distribution or invented tail extrapolation.

For F75L, expect a reference control fraction 0.234544 and a local supply slope
1.259147 V/V, compared with the saved backward native slope 1.258821 V/V.
Reference gate control contributes 1.225727 V/V of that local estimate,
while direct PMOS output conductance contributes 0.033419 V/V. The local
0.97-V extrapolation gives 0.462272 V versus the already observed 0.462088 V.

At 0.97 V, expect F75S/F75L means −38.456/−39.504 mV, sample SDs
32.252/24.066 mV, RMS errors 50.170/46.245 mV and MAEs 41.685/40.403 mV.
The fixed ±25-mV center counts are 165/147 of 521. The other four original
gain/bandwidth/current checks pass for these two circuits in both supply
conditions, so their five-check counts equal their center counts here.
The A75 0.97-V cohort retains one gain/signal-band failure as well as its
center failures; the analysis does not silently generalize the F75 equality.

## Check a predeclared capacitor choice

The second part supports chapters 5 and 8. It uses the same initial six-unit
D1 core as the previous learning-day package, with two later capacitor
conditions, **12 and 15 pF**. The [frozen plan](load-budget-plan.json) and
training hashes precede both new targets. The 500-ns deadline is explicitly
illustrative; it is separate from the unchanged finite DC tracking error.

The primary forecasts used the larger known 5/25-pF upper-bracket time per
farad. They predicted 12 pF would meet and 15 pF would miss the deadline in
both directions. The new observations preserve that result: 12-pF rise/fall
426.2–428.2 / 420.4–422.4 ns; 15-pF rise/fall 531.8–533.8 / 525.8–527.8 ns.
Those brackets are relative to the original input 50% crossings, within the
unchanged 1%-of-actual-OP-step band and original plateau horizons.

[Selected input metadata](data/interface-inputs.json) records both exact
circuits, physical units, original attempt identities, target times, recipes,
array columns and real parameter readbacks. The two NPZ files contain all
9,392 complete selected rows. There were two research attempts, no failed
acquisition, interruption or retry; the 15-pF deadline miss is a retained
physical requirement failure. Both target conditions are now exposed.

The five functions in `measurement_v1.py` are selected **verbatim** from the
earlier `learning.initial-d1-interface.measure.v1` implementation. Their
individual and module hashes are in the plan. `load_budget.py` adds only the
declared deadline classification; a spanning bracket stays unresolved, and
nonsettling or censored evidence cannot pass. It does not select a new load.

```bash
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B example/load_budget.py --output "$PWD/load-budget-001"
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B example/plot_budget.py --output "$PWD/load-budget-figures-001"
```

These commands need no original ZIP or SPICE installation: the exact two
new arrays, circuits and training-summary dependency are selected here.
They produce [the load-budget summary](load-budget-summary.json) and
`initial-d1-load-budget.svg/png`. Compare all three output DC values with the
known 5-pF endpoints: the observed differences are zero in stored float64
data. Bandwidth is 1.766517/1.413251 MHz. The 13.93-pF training estimate is
not a measured continuous maximum or a manufacturing tolerance.

## One advertised native 12-pF replay

Retrieve public APM v5.0.0, verify its exact commit and install ngspice **47**
separately. SPICE acquisition stays outside the static site build.

```bash
git clone --branch v5.0.0 --depth 1 https://github.com/ds54e/analog-process-models.git apm-v5
git -C apm-v5 rev-parse HEAD
ngspice --version
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B example/replay_budget.py --apm "$PWD/apm-v5" --binary /usr/local/bin/ngspice --output "$PWD/replay-12p-001"
```

The APM commit must be
`381517fda5107fabf98af7801d5a5103f38e230c`. Adjust only the binary path to the
installed ngspice 47. The script verifies the model hashes and tool version,
uses the exact selected 12-pF circuit and six saved zero-raw units, and writes
a new relocated realization, native logs and `result.json`. It neither imports
an ADS adapter nor accesses a private repository or legacy runtime.

The public API uses explicit terminal-voltage/current recipes with values
already present in the exact SPICE body. MOS model/W/L/m/nf/DELVTO/MULU0 are
read back by APM before, during and after the analyses. The three resistor
values are checked from observed DC voltage/current ratios; Cload is checked
from complex AC current/voltage and integrated transient charge. Expected
parameter values are not inserted as fake native readbacks.

The replay compares the three OP outputs, 361 AC points and the complete
adaptive transient with the saved 12-pF observation. It also checks the PWL,
maximum timestep, terminal KCL/domain, independent-endpoint settling bracket
and deadline result. This one nominal replay does not qualify a statistical
cohort, the entire capacitance interval or a real pass-gate interface.

## Checks and remaining limits

The [evidence record](evidence.json) retains the actual clean-directory check:
six summary/CSV outputs and three SVGs were regenerated byte-for-byte from
the selected ZIP and anonymously retrieved public inputs. The advertised
12-pF replay matched all three saved OP outputs, 361 AC output values and
4,333 transient output values with zero observed difference in this run.
The check used a fresh Python environment and public APM clone with the
existing ngspice 47 on the execution host; it was not a new OS installation
or an independent scientific replication.

The nodal matrix retains source degeneration, body effect at the NMOS bank,
finite Rf/input loading, PMOS output conductance and the finite reference
resistor. Its independently assembled linear equations agree with the
reduced expressions to floating-point precision. The original backward supply
steps and finite diagnostic perturbations are observations of the native
model; the equation omits leakage derivatives and capacitance. Their largest
local discrepancy across the selected input/current/gate-drop/supply paths
is below 0.03%. This is useful local agreement, not a continuous supply-domain
or frequency-domain guarantee.

Each output error is recalculated from the fixed 0.5-V target. The code checks
all original assignments, numerical/domain flags and five frozen performance
tests, reproduces the old mean/SD/RMS/MAE/count summaries, and separately checks
`RMS² = mean² + (n−1)/n × sample_SD²`. All valid physical failures are kept.
The mean and shape both change between F75S and F75L; the empirical comparison
does not isolate a variance-only intervention.

The supply/loss part acquires no new observations; the capacitor part adds
exactly its two predeclared nominal conditions. Neither part includes
per-sample trim, recentering, sample extension, a relaxed threshold, new
passive statistics, self-bias diagnosis or an LDO result.
Internal calculation/render checks are not prior human
review or independent scientific replication. Read the circuit choices and
costs in the [English route](https://ds54e.github.io/analog-design-notes/learn/index.html), with a
[Japanese overview](https://ds54e.github.io/analog-design-notes/ja/learning-path.html).
