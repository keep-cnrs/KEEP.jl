# **Web version of the jtcam-tex paper: plan**

Branch `jtcam-web` (from `jtcam-tex`). Everything new lives in `applications/jtcam-tex/web/`.

## Principle

Two sources only, never duplicated:

- the text: `main.tex` and `sections/*.tex` (the submission source);
- the results: `main.jl` and `src/*.jl` (every figure, table and quoted number).

The web site is a build product of both. It is never edited by hand.

## Decisions (jb, 2026-11-28)

1. One page per section (not one long page).
2. The web build reruns the full `main.jl` (about 15 min in CI): what the site shows is what the code it shows produces.
3. Hosting: subfolder `jtcam/` of KEEP's `gh-pages` (KEEP's own docs untouched).
4. `\alter` notes are removed from the web version.
5. Code granularity: a collapsible "code that produced this" block under each figure and table, plus full code pages.

## Pipeline

1. **Results.** Run `main.jl` (all sections, not `--fast`), as the `figures` job does. Then compare the regenerated `tables/` with the committed ones and fail on a difference (host-dependent values, such as `\SigmaGap`, get an explicit tolerance).
2. **Text: tex to md.** Pandoc (from `pandoc_jll`, so the same binary locally and in CI) reads `main.tex`, follows `\input`, and expands the simple `\newcommand` macros, including the generated numbers of `tables/numbers*.tex`. A prelude of web-only macro definitions (`web/macros.tex`) defines what pandoc cannot load from the `.sty` files (`\alter` as empty, `\vq`, `\dq`, ...). Lua filters (`web/filters/`) handle the rest:
   - `\SI`, `\ang` (siunitx) inside math: rewritten as `\text{}` units (MathJax has no siunitx);
   - cross-references: pandoc emits links with `reference-type` attributes; the filter numbers sections, figures, tables and equations as the PDF does, across pages;
   - one page per section: split at the level-1 headers;
   - figures: point to the generated PNGs; the two TikZ figures (`figs/cone_eight.tex`, `figs/configuration.tex`) are compiled to SVG in CI;
   - code markers (below).
3. **Maths.** MathJax 3 in Documenter, with the `physics` extension (`\vb`, `\pdv`, `\qty`, ... render unchanged) and AMS equation tags.
4. **Code shown with the results.**
   - A TeX comment next to each figure or table, `% web-code: src/analysis.jl#convergence`, names a region of the Julia source delimited by `#region convergence` ... `#endregion`. The filter inserts that code, read from the file at build time, as a collapsible block under the figure. The comments have no effect on the PDF.
   - Code pages: `main.jl` and `src/*.jl` through Literate.jl, plus `Project.toml`, `Manifest.toml`, the Julia version, the commit hash and the command to rerun.
   - API: the KEEP functions the paper calls (`LimitCycle.compute_limit_cycle`, `Optimization.optimize`, `SteadyState`, `PointMass4.dynamics`) through `@docs`, linked to the source at the paper's commit.
5. **Citations.** DocumenterCitations.jl with the same `refs.bib`; the filter turns pandoc's `[@key]` into `[key](@cite)`.
6. **CI.** A `web` job in `.github/workflows/jtcam-tex.yml`, on `workflow_dispatch` and tags: results, check, pandoc, `makedocs`, `deploydocs(dirname = "jtcam")`. Pull requests get a build without deployment.

## Phases

- **P1** (started): pandoc prototype, inventory of what breaks.
- **P2**: Documenter skeleton (`web/make.jl`, `web/Project.toml`), MathJax physics, citations, one page per section.
- **P3**: code regions and `% web-code:` markers, collapsible blocks.
- **P4**: CI job and deployment to `gh-pages/jtcam/`.
- **P5**: the check that regenerated results equal the committed ones.

## P1 findings (pandoc 3.14 from pandoc_jll, on 8344b6b)

- `pandoc main.tex -f latex -t markdown` runs without error and produces about 2150 lines.
- Handled by pandoc as is: `\input` of the sections, the generated number macros (e.g. the cycle periods appear as values), `\ac{}` (printed as the short form), `\cite` (as `[@key]`), `\ref`/`\eqref` (as links with `reference-type`), figures (as `<figure id=...>`), section labels (as header ids).
- To handle:
  - `\SI{...}{...}` inside math is left as is (13 occurrences): siunitx filter or MathJax macro.
  - `\alter` expands to garbage (`latex@warningalter [...]`): define it as empty in the web prelude.
  - 8 warnings for `\vb`/`\vu` in figure captions (alt text only; MathJax renders the caption itself).
  - Figures: only some `\includegraphics` came out as images; check the PDF-only ones and point all of them to the PNGs.
  - `\ac{}`: first use should give the long form (glossaries behaviour).
