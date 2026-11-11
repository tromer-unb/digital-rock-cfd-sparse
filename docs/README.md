# Documentation

Suggested reading order:

| # | document | for whom | what it contains |
|---|---|---|---|
| 1 | [../README.md](../README.md) | everyone | what is in the repository, quick start, reproduction table |
| 2 | [TUTORIAL.md](TUTORIAL.md) | first-time users | step by step: install, first run on an image of the paper, read the outputs, reproduce the paper, use your own images |
| 3 | [PAPER_EXPLAINED.md](PAPER_EXPLAINED.md) | readers of the paper | the method and the results in plain language, with pointers to code and data |
| 4 | [PROVENANCE.md](PROVENANCE.md) | reviewers, auditors | every figure, table and quoted number → result file → script |
| 5 | [REPRODUCIBILITY.md](REPRODUCIBILITY.md) | anyone re-running | environment, expected results, tolerances, resources |
| 6 | [AUDIT.md](AUDIT.md) | auditors | rules followed, checklist, worked example |

Documentation of the parts of the repository:

| folder | README |
|---|---|
| `PercoFlow/` | [PercoFlow/README.md](../PercoFlow/README.md) and [PercoFlow/docs/](../PercoFlow/docs/README.md): theory, inputs and outputs, OpenFOAM, benchmarks, glossary, code map, tutorial on synthetic examples |
| `data/` | [data/README.md](../data/README.md): every folder of results, naming, what is not stored |
| `reproduce/` | [reproduce/README.md](../reproduce/README.md) and [scripts_as_run/README.md](../reproduce/scripts_as_run/README.md) |
| `paper/` | [paper/README.md](../paper/README.md): figure files ↔ figure numbers |
| `legacy/` | [legacy/README.md](../legacy/README.md): the first release of this repository |
