# Analog Design Notes

Read the [site](https://ds54e.github.io/analog-design-notes/), start the
[eight-chapter learning route](https://ds54e.github.io/analog-design-notes/learn/),
or use the [Japanese reading overview](https://ds54e.github.io/analog-design-notes/ja/learning-path.html).

Learning edition **v1.1** connects MOS operating points, common-source stages,
mirrors, differential pairs/OTAs, feedback, errors and supply/startup to the
remaining interfaces of an LDO. Its [original selected evidence package](examples/learning-path-v1/README.md)
contains four archived nominal examples, current-balance calculations and eight
regenerable graphs. The [learning-day supplement](examples/learning-day-v1/README.md)
adds worked sizing/mirror data and a separately planned same-initial-core
5/25-pF receiving-port comparison with five new graphs and a native replay.
These are learning increments; the five
studies below retain their original histories and exact downloads.

The power-sequencing study explains input delivery/absorption and supply
charging peaks in 59 nominal static/history/control observations. It preserves
the driven-clamp boundary and finite-history limits. See the
[guide](https://ds54e.github.io/analog-design-notes/studies/power-sequencing/)
and [four small native examples](examples/power-sequencing/README.md).

The finite-feedback study connects lower active-load error to input-current,
resistor, noise and dynamic costs. Its frozen 521-coordinate confirmation and
24 preselected larger-signal cases preserve the limits of the improvement.
See the [guide](https://ds54e.github.io/analog-design-notes/studies/finite-feedback/)
and [public example](examples/finite-feedback/README.md).

The passive-sensitivity study asks when absolute or ratio resistor error changes
the benefit of added current. Its 5,568 new conditions reuse the known MOS cohort;
the passive scenarios are synthetic, not foundry statistics or fresh population
confirmation. See its [guide](https://ds54e.github.io/analog-design-notes/studies/passive-sensitivity/)
and [public example](examples/passive-sensitivity/README.md).

The first study compares gain, output operating point and untrimmed accuracy
in four finite amplifiers. It includes a frozen 464-coordinate paired study,
3,712 condition observations, exact circuit bodies, compact data and saved
Local physical parameters. Model limitations and valid failures are retained.

Start with the [example README](examples/gain-operating-point-accuracy/README.md)
to regenerate graphs or replay selected circuits. Native replay uses separately
retrieved APM v5.0.0 and ngspice 47. The static site build never runs SPICE.

Original prose/plots/data: CC BY 4.0. Original code/circuit examples: MIT.
See [LICENSE.md](LICENSE.md) and [third-party notices](THIRD_PARTY_NOTICES.md).
These releases have internal AI checks, without prior human or peer review.

The editable canonical manuscripts are maintained in the private research
repository. This repository receives an explicit selected publication set.
Reader corrections belong in Issues; they are applied to that canonical source
before export. The immutable research tags are
`gain-operating-point-accuracy-v1.0`, `passive-sensitivity-v1.0`,
`active-allocation-v1.0`, `finite-feedback-v1.0` and `power-sequencing-v1.0`.
The original edition tag `learning-path-v1.0` is retained; the new edition tag
is `learning-path-v1.1`. Versioned
downloads are static checked artifacts; the first release's ZIP bytes are preserved.

Maintainers render with Quarto **1.10.18**:

```bash
python3 tools/build_site.py
quarto render site
python3 tools/check_site.py --site site/_site
```

The `Publish Pages` workflow is dispatched on `main` after a selected release
is ready. It builds only this public repository, records its exact commit in
every page and `revision.json`, and deploys the Pages artifact. No private
repository, model acquisition, credentials file or research runtime is bundled.

The active-area allocation guide adds a fixed 598-coordinate comparison of
four same-area active redesigns with two matched-current baselines. It preserves
failed nominal design attempts and inconclusive minimum-effect decisions.
See the versioned [guide](https://ds54e.github.io/analog-design-notes/studies/active-allocation/)
and its [reproduction example](https://ds54e.github.io/analog-design-notes/studies/active-allocation/reproduce.html).
