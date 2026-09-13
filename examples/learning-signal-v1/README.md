# Reproduce the local-gain and finite-signal explanation

This learning-signal-v1.0 supplement accompanies learning edition v1.3.
It connects chapters 2 and 5: the derivative
of a complete circuit, its finite-signal harmonics, internal gate/source
motion and the current its input source must supply or absorb.

Every scientific observation here was already acquired in the five-study
work. The reanalysis draws no new physical sample, changes no original
confirmation or threshold, and retains R75's valid 20-mV distortion failure.
The additional selected channels are saved internal-node observations from
those same attempts. Their publication does not make them fresh confirmation.

Download [the selected package](https://ds54e.github.io/analog-design-notes/downloads/learning-signal-v1.0.zip),
or inspect [the public files](https://github.com/ds54e/analog-design-notes/tree/main/examples/learning-signal-v1).
The original studies and all earlier versioned downloads remain unchanged.

## Exact public inputs and two different measurements

The external input is the preserved
[feedback-study-v1.0.zip](https://ds54e.github.io/analog-design-notes/downloads/feedback-study-v1.0.zip), SHA256
`348d3f2c8cb3ef41b9c968f4bbb68c011c40cfb64ddda4138d2ff115c5e95bc6`.
Its original tag is finite-feedback-v1.0 at public commit
`721a59b9c33146b6ff1d8662e6b9a5076aafe448`.

[The source record](data/source.json) lists eleven exact members, their
hashes, six original nominal sine identities, model hashes, circuit
parameters and zero-raw physical units. R75 and F75L each contribute
10/20-mV-peak cycles at 10 kHz and 1 pF; the two original 20-mV numerical
refinements remain separate cases. All use 1 V, 26.85 °C, 0.45-V nominal input,
100-kΩ output load and the original finite source/bias networks.

The original public table has one final cycle on a 1,001-point uniform grid,
0.9–1.0 ms including the repeated endpoint. The new port-cycle array selects
the actual gate/source voltages and NMOS terminal-current sums from the same
native traces. All six already public columns match their original
interpolated values exactly. The remaining full traces stay with the
preserved research evidence; no fresh transient was needed to select them.

The new analysis uses the **first 1,000 points**, excluding the duplicate
endpoint, and resolves peak complex harmonics 1–10. It compares the FFT with
an independently assembled sine/cosine least-squares basis. The original
requirement used **eight cycles** and its own stationarity test. Those
original metrics and the new one-cycle decomposition are separate fields;
one saved cycle cannot independently prove stationarity.

[The static source record](data/static-source.json) and
[eight selected DC points](data/dc-points.csv) retain R75's ±0.2/±0.1-mV
input perturbations, F75L's ±0.1-mV perturbations and their two original
nominal centers. Models, physical units, input settings, native OP values
and exact circuit bytes were checked. R75's older perturbation bodies retain
their 10-mV sine text; the OP recipes explicitly set DC input values.
The later 5-mV nominal center has the same MOS inventory/raw values and
passives. The different body, recipe and bound-realization identities are
preserved, not replaced by a common identifier.

## Recalculate the explanation and graphs

Use a new directory with the extracted example, Python 3.9 or later, and
the supplied pinned requirements. These calculations need no SPICE or APM
installation. Retrieve the original ZIP from the public URL above.

~~~bash
python3 -m venv .venv
.venv/bin/python -m pip install -r example/requirements.txt
curl --fail --location https://ds54e.github.io/analog-design-notes/downloads/feedback-study-v1.0.zip --output feedback-study-v1.0.zip
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B example/analyze.py --feedback-zip feedback-study-v1.0.zip --output "$PWD/analysis-001"
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B example/plot.py --feedback-zip feedback-study-v1.0.zip --output "$PWD/figures-001"
~~~

Existing output destinations are refused. The analysis writes five CSV
tables and a summary JSON: all harmonics, cycle metrics, two retrospective
extrapolations, the local resistor/MOS network, and static finite differences.
The plots are local-gain-and-curvature and feedback-gate-motion, in SVG/PNG.

Expect static second derivatives of approximately −50.1781 V⁻¹ for R75
and −7.57495 V⁻¹ for F75L. The R75 two-step change is about 0.00047%;
F75L has one saved static step. The corresponding quadratic 20-mV THD
estimates are 4.182%/0.631%, versus original observations 4.156%/0.671%.
An additional explicitly retrospective 10-to-20-mV sine extrapolation is
retained; it is not a forecast evaluated on an unseen target.

At 20 mV, expect gate/source/gate-minus-source fundamental peaks
19.998/3.295/16.703 mV in R75 and 7.290/1.816/5.474 mV in F75L.
The voltage differences are taken before their Fourier decomposition.
The actual gate voltage independently checks the measured-input-current
resistor identity; source voltage checks the summed native source current.

## One native replay of the known F75L example

Retrieve public APM v5.0.0 and install ngspice **47** separately. The APM
commit must be `381517fda5107fabf98af7801d5a5103f38e230c`. SPICE acquisition
is outside the static site build.

~~~bash
git clone --branch v5.0.0 --depth 1 https://github.com/ds54e/analog-process-models.git apm-v5
git -C apm-v5 rev-parse HEAD
ngspice --version
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B example/replay.py --apm "$PWD/apm-v5" --feedback-zip feedback-study-v1.0.zip --binary /usr/local/bin/ngspice --output "$PWD/replay-F75L-001"
~~~

Adjust only the binary path to the installed ngspice 47. The script uses
[the exact 20-mV circuit](circuits/F75L-20m.cir) and its five saved zero-raw
units through the public APM recipe API. It records a new relocated
realization and native run; it does not import an ADS adapter or access a
private repository/runtime. The [finite replay plan](replay-plan.json) fixes
the three known DC inputs, central AC sweep, original 1-ms transient and
comparison tolerances before this replay.

The native output checks the old static derivative and complete selected
cycle. The original eight-cycle function is selected verbatim, with its
module/function identities in the plan. Actual model/W/L/m/nf/DELVTO/MULU0
are read back before, during and after the analyses. Six resistor values
are inferred from measured voltage/current, including feedback current from
gate KCL. The capacitor is checked through AC current/voltage and integrated
charge. Expected parameter values are not inserted as native readbacks.

This is a reproduction of known nominal conditions. It does not acquire a
fresh confirmation or extend the demonstrated amplitude range.

## Checked public reproduction

On 2026-09-13, a new directory and fresh Python environment used the selected
package, anonymously retrieved original feedback ZIP and public pinned APM.
All six summary/CSV files and both SVG graphs regenerated byte-for-byte.
The native F75L replay passed on the existing ngspice 47 installation:
three OP records, 401 AC points and 10,011 transient points. The seven selected
cycle voltage/current columns and original eight-cycle THD matched the saved
values exactly; the DC curvature differed by 0.0031%, within its predeclared
0.05% comparison tolerance. Actual MOS/passive readback and KCL/charge checks
also passed. The original measurement and its thresholds were unchanged.

This is same-host reproduction with fresh dependencies and public inputs,
not a clean operating-system installation or independent scientific replication.
Detailed numerical differences and input identities are in the evidence record.

## Interpretation and limits

The static derivative belongs to the whole circuit after its current
balances are solved. The apparent curvature inferred from a finite sine
is a different estimate. For F75L it changes with amplitude: the 20-mV
inferred magnitude is about 5.5% above the local DC difference. Higher
orders and memory remain visible, rather than being fitted away.

The quadratic distortion allowance is a first screening estimate.
It does not establish a continuous input range, a rail-to-rail interface,
performance of a new physical population, or a regulator. The algebraic
curvature multiplier uses a specified gate-driven diagnostic fixture,
including Rf and source degeneration; it does not isolate Rf experimentally
or measure a return ratio/phase margin.

F75L's input source must deliver and absorb current. Its smaller distortion
coexists with higher noise, more internal resistance and lower bandwidth
than R75 under the original same-function comparison. The time-average
shift of one driven waveform is separate from the original DC output-center
distribution across saved Local realizations.

The [evidence record](evidence.json) distinguishes the saved-data checks and
the public-input reproduction checks. These are internal checks on the
execution host, without prior human review or independent replication.
The original study's separate
[native reproduction](https://ds54e.github.io/analog-design-notes/studies/finite-feedback/reproduce.html) remains available
under its unchanged inputs and definitions.
