# Reproduce the initial-OTA DC receiving-load example

This example changes the stationary current demanded by initial D1's receiver.
It keeps the same six zero-raw physical units, 0.3-V source, finite bias/source
network and 5-pF capacitor. Rload changes from the known 1 MΩ to 100 kΩ.
Chapters [4](https://ds54e.github.io/analog-design-notes/learn/04-differential-pair-and-ota.html#a-dc-receiver-needs-current-after-charging-has-ended),
[5](https://ds54e.github.io/analog-design-notes/learn/05-feedback-and-response.html#predict-a-finite-dc-receiving-current-change)
and [8](https://ds54e.github.io/analog-design-notes/learn/08-toward-an-ldo.html) explain the resulting current,
precision and response tradeoff.

Download [the selected package](https://ds54e.github.io/analog-design-notes/downloads/learning-dc-load-v1.0.zip)
or inspect [its public files](https://github.com/ds54e/analog-design-notes/tree/main/examples/learning-dc-load-v1).
This package has one usable new OP/AC condition. It also documents the earlier
wrong-condition acquisition and the correction; its completion marker was
never accepted as evidence for the intended load.

## Inputs, versions and exposed conditions

The [exact circuit](circuits/D1-100k-load.cir) uses VDD = 1 V, 26.85 °C,
Rsource = 10 kΩ, Rbias = 273692.63541213644 Ω, Cload = 5 pF and Rload =
100 kΩ. Every MOS drain/gate/source/body terminal and diagnostic sensor is
explicit. Geometry, device inventory and saved raw coordinate pairs are the
same as the initial core. Unused Vin PWL text is removed for this stationary
fixture; Vin is DC 0.3 V, AC 1. A unit AC magnitude is a local normalization,
not a claimed valid 1-V large-signal excursion.

The [corrected machine plan](plan.json) freezes one OP and 361 AC points,
60 per decade from 1 kHz to 1 GHz. It includes every source/element override,
physical unit, model hash and numerical tolerance. The [observation index](data/index.json)
preserves the actual start, request/realization identity, native deck and
output hashes. Complete arrays and actual MOS/passive readbacks were checked
before measurement. No new noise, transient, population or load search was run.

The prior package is [learning-day-v1.0.zip](https://ds54e.github.io/analog-design-notes/downloads/learning-day-v1.0.zip),
SHA256 `ddfbad7a06c9c45a94fa4f0379a421816f04d3fa1c56bd21ec8a09ac99f13b15`.
[source.json](data/source.json) identifies three exact members: the initial
configuration, D1-5p circuit and D1-5p observations. Only its central OP/AC
are used here. Its full trajectory and the separate 25-pF cases retain their
original definitions and exposure history. The archived selected-D1 design
is not substituted for this initial six-unit core.

### Retained construction error

The [original input plan](original-input-plan.json) changed the circuit body
and expected fields to 0.3 V / 100 kΩ, but accidentally retained the earlier
per-analysis overrides Vin = 0.7 V and Rload = 1 MΩ. The simulator applied
those overrides. The unchanged passive-parameter check rejected the result;
actual source voltage also contradicted the intended condition.

The [correction record](input-correction.json) retains its identities and
actual readback. Every selected OP/AC value repeated the already exposed
high-command observation in the [common-mode package](https://ds54e.github.io/analog-design-notes/studies/learning-common-mode-v1/reproduce.html).
It was a numerically complete acquisition at the wrong condition. It provided
no observation of the intended 100-kΩ target.

A separately versioned freeze corrected only those two override values in
both analyses. The body, six physical units, original prior forecast,
measurement code and tolerances remained byte-identical. Preparation now
checks agreement among body, expected fields and all overrides. Exactly one
corrected acquisition followed; it passed. There were **two acquisitions in
this item, one usable intended condition and no further scientific retry**.
Original attempts remain preserved; the first completion marker is not
relabeled as valid 100-kΩ evidence.

## Recompute the saved observations

Extract the selected ZIP into a new directory. It contains `example/`, with
the exact scripts, compact arrays, source hashes and retained license notices.
Use Python 3.9 or later and the established pinned requirements:

~~~bash
python3 -m venv .venv
.venv/bin/python -m pip install -r example/requirements.txt
curl --fail --location https://ds54e.github.io/analog-design-notes/downloads/learning-day-v1.0.zip --output learning-day-v1.0.zip
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B example/prior.py --source-zip learning-day-v1.0.zip --output "$PWD/prior-001"
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B example/analyze.py --source-zip learning-day-v1.0.zip --output "$PWD/analysis-001"
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B example/explain.py --source-zip learning-day-v1.0.zip --output "$PWD/network-001"
~~~

Expected files are `prior.json`, `summary.json`, `operating-point.csv`,
`network.json` and `network.csv`. Existing output directories are refused.
The two CSV writers explicitly use LF record separators. The [manifest](data/manifest.json)
binds selected data; measurements also check their frozen script and inputs.
No new plot is needed for this single-condition table comparison; the reading
route retains its existing graph sources and reproduction packages.

Actual resistor parameters and independent DC V/I agree, including the
100-kΩ receiving load. All MOS and node current balances, source voltages
and terminal spans pass the retained checks. The capacitor is independently
checked over 100 kHz–10 MHz using
`Icap = -i(VtM2d)-i(VtM4d)-Vout/Rload+i(Vfeedback)` and
`C = Icap/(j*2*pi*f*Vout)`. Its DC current is zero; its AC current is observed.
Bandwidth retains the first downward −3.000-dB crossing relative to the
observed 1-kHz gain, interpolated in dB/log frequency. It is not the five-study
half-power interpolation definition or a measured return ratio.

## Keep the prior forecast separate from its explanation

The prior calculation uses the known 1-MΩ operating point. Two node-equation
constructions give a local output-current transimpedance of 7490.861 Ω.
For added conductance 9 µS, its affine prediction gives Vout = 0.279606 V,
real 1-kHz gain = 0.921211 and bandwidth = 4.527397 MHz. The new observations
are 0.279412 V, 0.918230 and 4.401678 MHz. The −0.193-mV output discrepancy
and −2.78% bandwidth discrepancy are retained.

The separate `explain.py` calculation is retrospective. With the newly
observed gm/gmb/gds it gives gain 0.918337 and a −3.000-dB bandwidth estimate
of 4.399545 MHz from the external capacitor alone. These do not replace the
saved prior prediction. Body derivatives are retained; leakage derivatives
and intrinsic capacitances are omitted in the local network.

A separate stationary terminal-power sum includes actual gate/body currents,
all three resistors and signed input/supply power. Its residual is below
1 fW. Load current increases 0.298→2.794 µA while supply current changes
16.011→16.029 µA: increased P2 sourcing and reduced mirrored N4 sinking
provide the receiver current together. Load power rises 0.089→0.781 µW,
with changed internal dissipation and lower output accuracy. This is not a
transient-energy, efficiency-optimum or integrated-regulator claim.

## Replay the known corrected condition natively

The known-input [replay plan](replay-plan.json) was saved after the scientific
condition was exposed. It is reproduction, not fresh confirmation. Obtain the
public APM source at its pinned commit and ngspice 47:

~~~bash
git clone --branch v5.0.0 --depth 1 https://github.com/ds54e/analog-process-models.git apm-v5
git -C apm-v5 rev-parse HEAD
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python -B example/replay.py --source-zip learning-day-v1.0.zip --apm "$PWD/apm-v5" --binary /usr/local/bin/ngspice --output "$PWD/native-replay-001"
~~~

The reported source commit must be
`381517fda5107fabf98af7801d5a5103f38e230c`.
The script verifies body/recipe agreement before native execution, the pinned
model bytes and ngspice version. It relocates the same six saved zero-raw units
and records both realization identities. The public API replay uses terminal
voltage/current vectors, so independent resistor V/I and AC C readback support
its passive checks. The scientific saved observations retain their additional
native parameter/derivative vectors. All selected terminal observations are
compared at the original OP/AC axes using the declared replay tolerances.

The clean-directory check on 2026-09-13 retrieved the old ZIP and pinned APM
source anonymously, installed the declared Python requirements and replayed
this known condition. All selected terminal voltages/currents matched the
saved observations exactly at one OP and 361 AC points. The prior, measurement
and retrospective calculation regenerated all five JSON/CSV result files with
exact bytes. The [evidence record](evidence.json) retains input and replay
identities. This same-host check is not independent scientific replication.

This one diagnostic resistor and capacitor do not model a complete receiving
MOS gate, nonlinear charge, overload recovery or LDO. The closed response alone
does not establish loop margins. No source-profile recalibration, device
redraw, resistor retuning or claim of prior human review is involved.
