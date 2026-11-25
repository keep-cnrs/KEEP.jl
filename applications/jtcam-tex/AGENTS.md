# AGENTS.md — Rocking Arm Model for Kite Power Generation

## Build

```sh
latexmk -pdf -outdir=latex_out _main.tex
```

Class option `OptSoft` requires **biber** + `software-biblatex`.  
`makeglossaries` needed for acronyms (`myacronyms.sty`).

## Structure

- `_main.tex` — entry point. Only `sections/section1-modelling.tex` included; appendices (`sectionA-10D.tex`, `sectionB-explicit-ddq.tex`) commented out.
- `sectionAA-system-gemini.tex` — standalone TikZ figure doc, compiles independently.
- `jtcam_preprint.cls` / `simple_jtcam_preprint.cls` — journal templates.
- `latex_out/` — build output.
- Two `.bib` files: both contain placeholder entries (Salam, Pontzen, Kutz) — replace before submission.

## Physics package notation

Already loaded (`_main.tex:6`). Prefer these over raw `\frac`:
- `\vb`, `\vu` — bold/unit vectors
- `\pdv{f}{x}`, `\pdv*{f}{x}` — partial derivatives (star = displayed)
- `\dv{f}{t}`, `\dv*{f}{t}` — total time derivatives
- `\norm{...}`, `\abs{...}` — norms and absolute values
- Only use summations over components if 3rd order tensors are involved, otherwise use linear algebra. Also explore the possibility of using directional derivatives to avoid summations.

## Custom packages

- `mymacros.sty` — math shorthands: `\vq`, `\dq`, `\ddq`, `\R`, `\ahat`, `\rhat`, `\lhat`, `\OA`, `\AK`, `\OK`
- `myacronyms.sty` — ODE, DAE, AD, DoF

## Conventions

- `\alter{...}` — WIP note, renders brown. Used for TODOs and placeholders.
