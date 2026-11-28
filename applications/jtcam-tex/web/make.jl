# Web version of the paper: main.tex + the outputs of main.jl -> Documenter site.
#
#   julia --project=applications/jtcam-tex/web applications/jtcam-tex/web/make.jl [--aux=main.aux] [--deploy]
#
# Steps (nothing here is edited by hand, see web/PLAN.md):
#   1. numbers of sections, equations, figures and tables: from main.aux (the LaTeX build is
#      run here unless --aux is given), so that the site and the PDF agree;
#   2. text: pandoc (pandoc_jll) reads main.tex with a prelude of web macros, runs citeproc
#      on refs.bib and the Lua filter web/filters/paper.lua, which writes one page per
#      section in web/src/;
#   3. code: the `% web-code: file#region` comments of the .tex files name `#region` blocks
#      of the Julia sources (or whole TikZ files), shown under the figures and tables; code
#      pages show main.jl and src/*.jl in full;
#   4. Documenter builds web/build; --deploy pushes it to the `jtcam/` folder of gh-pages.
#
# Figures and tables are taken from figs/ and tables/ as they are (the CI job regenerates
# them with main.jl first, and web/check_results.jl compares the tables with the committed ones).

using Documenter
using JSON
using pandoc_jll
using KEEP

const WEB = @__DIR__
const ROOT = dirname(WEB)                    # applications/jtcam-tex
const REPO = dirname(dirname(ROOT))          # the KEEP.jl repository
const SRC = joinpath(WEB, "src")             # generated pages (not committed)
const GITHUB = "github.com/keep-cnrs/KEEP.jl"
const RELROOT = "applications/jtcam-tex"

const ARGS_DEPLOY = "--deploy" in ARGS
const ARG_AUX = let a = findfirst(startswith("--aux="), ARGS)
    a === nothing ? nothing : abspath(ARGS[a][7:end])
end

commit() = get(ENV, "GITHUB_SHA", readchomp(`git -C $REPO rev-parse HEAD`))
blob_url(file, lines=nothing) = "https://$GITHUB/blob/$(commit())/$RELROOT/$file" *
                                 (lines === nothing ? "" : "#L$(first(lines))-L$(last(lines))")

# ------------------------------------------------------------------ 1. LaTeX build, main.aux

"Compile the paper (latexmk) in a temporary directory; return the directory."
function latex_build()
    dir = mktempdir()
    @info "web: compiling main.tex for main.aux (latexmk)"
    cd(ROOT) do
        run(`latexmk -pdf -interaction=nonstopmode -silent -outdir=$dir main.tex`)
    end
    return dir
end

"Labels of main.aux: label => (num, kind)."
function parse_aux(path)
    labels = Dict{String,Any}()
    for line in eachline(path)
        m = match(r"^\\newlabel\{([^}]+)\}\{\{([^}]*)\}", line)
        m === nothing && continue
        lab, num = m.captures
        endswith(lab, "@cref") && continue
        k = match(r"\{(Doc-Start|section|subsection|subsubsection|appendix|equation|figure|table)[^}]*\}\{[^}]*\}\}\s*$", line)
        kind = k === nothing ? "other" : k.captures[1]
        kind = kind in ("section", "subsection", "subsubsection", "appendix", "Doc-Start") ? "section" : kind
        labels[lab] = Dict("num" => num, "kind" => kind)
    end
    return labels
end

# ------------------------------------------------------------------ 3. code shown with the results

"Lines of `file` between `#region name` and `#endregion` (exclusive), and their range."
function region(file, name)
    lines = readlines(joinpath(ROOT, file))
    i = findfirst(l -> strip(l) == "#region $name", lines)
    i === nothing && error("web-code: no `#region $name` in $file")
    j = findnext(l -> strip(l) == "#endregion", lines, i + 1)
    j === nothing && error("web-code: `#region $name` of $file has no `#endregion`")
    body = lines[i+1:j-1]
    ind = minimum((length(l) - length(lstrip(l)) for l in body if !isempty(strip(l))); init=0)
    return join((l[min(ind, length(l))+1:end] for l in body), "\n"), (i + 1):(j - 1)
end

slug(file) = replace(file, r"[/.]" => "-")
code_page(file) = "code/$(slug(file)).md"
lang(file) = endswith(file, ".jl") ? "julia" : "latex"

"The .tex files in the order of main.tex, with the generated tables expanded in place."
function tex_lines()
    out = String[]
    function expand(file)
        for l in eachline(joinpath(ROOT, file))
            m = match(r"^\s*\\input\{([^}]+)\}", l)
            if m !== nothing && !startswith(strip(l), "%")
                f = m.captures[1]
                f = endswith(f, ".tex") ? f : f * ".tex"
                startswith(f, "tables/numbers") || expand(f)
            end
            push!(out, l)
        end
    end
    expand("main.tex")
    return out
end

"label => code items (from `% web-code:`) and label => TikZ figure (from \\includestandalone)."
function scan_tex()
    webcode = Dict{String,Any}()
    standalone = Dict{String,String}()
    last_label, pending_sa = nothing, nothing
    for l in tex_lines()
        m = match(r"\\includestandalone(?:\[[^\]]*\])?\{([^}]+)\}", l)
        m !== nothing && (pending_sa = m.captures[1])
        m = match(r"\\label\{([^}]+)\}", l)
        if m !== nothing && !startswith(strip(l), "%")
            last_label = m.captures[1]
            if pending_sa !== nothing
                standalone[last_label] = pending_sa
                pending_sa = nothing
            end
        end
        m = match(r"^\s*%\s*web-code:\s*(.*)$", l)
        m === nothing && continue
        last_label === nothing && error("web-code marker before any \\label: $l")
        items = get!(webcode, last_label, Any[])
        for spec in split(strip(m.captures[1]))
            file, name = occursin('#', spec) ? split(spec, '#'; limit=2) : (spec, "")
            code, lines = name == "" ? (rstrip(read(joinpath(ROOT, file), String)), nothing) : region(file, name)
            push!(items, Dict("file" => file, "region" => name, "code" => code, "lang" => lang(file),
                "url" => blob_url(file, lines), "page" => code_page(file)))
        end
    end
    return webcode, standalone
end

# ------------------------------------------------------------------ figures

"TikZ figure figs/<name>.tex -> src/figs/<name>.svg (pdflatex + pdftocairo), or .png with gs."
function tikz_figure(name)
    dir = mktempdir()
    cp(joinpath(ROOT, "$name.tex"), joinpath(dir, "fig.tex"))
    cd(dir) do
        run(pipeline(`pdflatex -interaction=nonstopmode fig.tex`; stdout=devnull))
    end
    out = joinpath(SRC, name)
    if Sys.which("pdftocairo") !== nothing
        run(`pdftocairo -svg $(joinpath(dir, "fig.pdf")) $out.svg`)
        return "$name.svg"
    elseif haskey(ENV, "JTCAM_GS")  # local fallback: ghostscript
        run(`$(ENV["JTCAM_GS"]) -q -dSAFER -dBATCH -dNOPAUSE -sDEVICE=pngalpha -r200 -sOutputFile=$out.png $(joinpath(dir, "fig.pdf"))`)
        return "$name.png"
    end
    @warn "web: no pdftocairo (nor JTCAM_GS): TikZ figure $name not converted"
    return "$name.svg"
end

# ------------------------------------------------------------------ 2. pandoc

"""
Copy of the paper sources for pandoc, with the web-only rewrites of TEX_FIXES (constructs
that pandoc cannot read, such as the siunitx `S` columns). The sources are not modified.
"""
function pandoc_tree()
    dir = mktempdir()
    for f in readdir(ROOT)
        (endswith(f, ".tex") || endswith(f, ".sty") || endswith(f, ".bib")) && cp(joinpath(ROOT, f), joinpath(dir, f))
    end
    for sub in ("sections", "tables")
        mkpath(joinpath(dir, sub))
        for f in readdir(joinpath(ROOT, sub))
            endswith(f, ".tex") || continue
            s = read(joinpath(ROOT, sub, f), String)
            for (re, sub_) in TEX_FIXES
                s = replace(s, re => sub_)
            end
            write(joinpath(dir, sub, f), s)
        end
    end
    cp(joinpath(WEB, "macros.tex"), joinpath(dir, "web-macros.tex"))
    return dir
end

# siunitx S columns -> r (pandoc drops the whole tabular otherwise)
const TEX_FIXES = [
    r"\\begin\{tabular\}\{([^}]*)\}" => s -> replace(s, r"(?<=\{)[^}]*(?=\}$)" => c -> replace(c, 'S' => 'r')),
]

"Run pandoc + citeproc + web/filters/paper.lua; web/macros.tex is read first, so its definitions win."
function run_pandoc(datafile)
    dir = pandoc_tree()
    cd(dir) do
        run(`$(pandoc()) web-macros.tex main.tex -f latex -t html --math-method=mathjax -o $(devnull)
             --citeproc --bibliography=refs.bib -M link-citations=true
             -M webdata=$datafile -M outdir=$SRC --lua-filter=$(joinpath(WEB, "filters", "paper.lua"))`)
    end
end

# ------------------------------------------------------------------ extra pages

"Title, authors and keywords of main.tex."
function front_matter()
    tex = read(joinpath(ROOT, "main.tex"), String)
    title = match(r"\\title\{([^}]*)\}", tex).captures[1]
    authors = [replace(m.captures[1], r"\\orcid\{[^}]*\}" => "") for m in eachmatch(r"\\author\[[^\]]*\]\{(.*)\}", tex)]
    kw = match(r"\\keywords\{([^}]*)\}", tex)
    return (; title, authors, keywords=kw === nothing ? "" : kw.captures[1])
end

function write_index(pages, fm)
    io = IOBuffer()
    println(io, "# $(fm.title)\n")
    println(io, "*", join(fm.authors, ", "), "*\n")
    isempty(fm.keywords) || println(io, "**Keywords:** ", fm.keywords, "\n")
    println(io, """
    This is the web version of the paper, built from the same sources as the [PDF](assets/main.pdf):
    the text of `main.tex` and the results of `main.jl` (every figure, table and number quoted
    in the text). Under each figure and table, *Code that produced this* shows the Julia code
    that computed it; the [code pages](code/main-jl.md) give the full sources and
    [Reproducing the results](reproduce.md) the commands to rerun them.

    Built from commit [`$(commit()[1:7])`](https://$GITHUB/tree/$(commit())/$RELROOT).

    ## Contents
    """)
    for p in pages
        println(io, "- [$(p["title"])]($(p["file"]))")
    end
    write(joinpath(SRC, "index.md"), String(take!(io)))
end

const CODE_FILES = ["main.jl", "src/helpers.jl", "src/modelling.jl", "src/numerics_10d.jl",
    "src/analysis.jl", "src/optimization.jl", "figs/configuration.tex", "figs/cone_eight.tex"]

function write_code_pages()
    mkpath(joinpath(SRC, "code"))
    for f in CODE_FILES
        src = rstrip(read(joinpath(ROOT, f), String))
        write(joinpath(SRC, code_page(f)), """
        # `$f`

        [View on GitHub]($(blob_url(f))) (commit `$(commit()[1:7])`).

        ```$(lang(f))
        $src
        ```
        """)
    end
end

function write_reproduce()
    proj = rstrip(read(joinpath(ROOT, "Project.toml"), String))
    write(joinpath(SRC, "reproduce.md"), """
    # Reproducing the results

    Every figure (`figs/*.png`, `figs/*.pdf`), table (`tables/table*.tex`) and number quoted in
    the text (`tables/numbers.tex`, `tables/numbers_analysis.tex`) is written by `main.jl`. From
    the root of the [KEEP.jl repository](https://$GITHUB), at commit `$(commit())`:

    ```sh
    git checkout $(commit())
    julia --project=$RELROOT -e 'import Pkg; Pkg.instantiate()'
    julia --project=$RELROOT $RELROOT/main.jl            # all sections, about 15 min
    julia --project=$RELROOT $RELROOT/main.jl analysis   # one section: modelling, 10d, analysis, optimization
    ```

    Julia $(VERSION) built this site. The package versions are pinned by
    [`Manifest.toml`]($(blob_url("Manifest.toml"))); the environment is

    ```toml
    $proj
    ```

    The paper itself: `latexmk -pdf main.tex` in `$RELROOT` (biber and `makeglossaries` needed).
    """)
end

const API = [
    "KEEP.PointMass4" => ["KEEP.PointMass4.dynamics", "KEEP.PointMass4._accelerations_tension",
        "KEEP.PointMass4.compute_line_tension", "KEEP.PointMass4.jacobian"],
    "KEEP.LimitCycle" => ["KEEP.LimitCycle.compute_limit_cycle", "KEEP.LimitCycle.refine_shooting",
        "KEEP.LimitCycle.all_limit_cycles", "KEEP.LimitCycle.average_power"],
    "KEEP.SteadyState" => ["KEEP.SteadyState.all_steady_states"],
    "KEEP.Optimization" => ["KEEP.Optimization.optimize"],
    "KEEP.PointMassPara" => ["KEEP.PointMassPara.optimal_torque_slope"],
]

function write_api()
    io = IOBuffer()
    println(io, "# KEEP functions used by the paper\n")
    println(io, "The docstrings of the main functions of [KEEP.jl](https://$GITHUB) called by `main.jl`.\n")
    for (mod, fs) in API
        println(io, "## `$mod`\n\n```@docs")
        foreach(f -> println(io, f), fs)
        println(io, "```\n")
    end
    write(joinpath(SRC, "api.md"), String(take!(io)))
end

# ------------------------------------------------------------------ main

function build()
    rm(SRC; force=true, recursive=true)
    mkpath(joinpath(SRC, "figs"))
    mkpath(joinpath(SRC, "assets"))

    texdir = ARG_AUX === nothing ? latex_build() : dirname(ARG_AUX)
    aux = ARG_AUX === nothing ? joinpath(texdir, "main.aux") : ARG_AUX
    labels = parse_aux(aux)
    pdf = joinpath(texdir, "main.pdf")
    isfile(pdf) && cp(pdf, joinpath(SRC, "assets", "main.pdf"))
    cp(joinpath(WEB, "assets", "paper.css"), joinpath(SRC, "assets", "paper.css"))

    for f in readdir(joinpath(ROOT, "figs"))
        endswith(f, ".png") && cp(joinpath(ROOT, "figs", f), joinpath(SRC, "figs", f))
    end
    webcode, standalone = scan_tex()
    for (lab, name) in standalone
        standalone[lab] = tikz_figure(name)
    end

    datafile = joinpath(mktempdir(), "webdata.json")
    write(datafile, JSON.json(Dict("labels" => labels, "webcode" => webcode, "standalone" => standalone)))
    @info "web: pandoc ($(length(labels)) labels, code for $(length(webcode)) figures and tables)"
    run_pandoc(datafile)

    pages = JSON.parsefile(joinpath(SRC, "pages.json"))
    rm(joinpath(SRC, "pages.json"))
    rm(joinpath(SRC, "front.html"))
    fm = front_matter()
    write_index(pages, fm)
    write_code_pages()
    write_reproduce()
    write_api()

    is_app(p) = p["appendix"]
    main_pages = [p["title"] => p["file"] for p in pages if !is_app(p) && p["file"] != "references.md"]
    app_pages = [p["title"] => p["file"] for p in pages if is_app(p)]
    has_refs = any(p -> p["file"] == "references.md", pages)

    mathengine = Documenter.MathJax3(Dict(
            :loader => Dict("load" => ["[tex]/physics", "[tex]/mathtools"]),
            :tex => Dict(
                "inlineMath" => [["\\(", "\\)"]],
                "displayMath" => [["\\[", "\\]"]],
                "processEnvironments" => true,
                "tags" => "none",
                "packages" => Dict("[+]" => ["physics", "mathtools"]),
            ),
        ), true)

    makedocs(;
        root=WEB, source="src", build="build",
        sitename="Rocking arm model for kite power generation",
        modules=[KEEP], checkdocs=:none, warnonly=true, remotes=nothing,
        format=Documenter.HTML(;
            prettyurls=false, mathengine, size_threshold=nothing, size_threshold_warn=nothing,
            repolink="https://$GITHUB", assets=["assets/paper.css"],
            edit_link=nothing, inventory_version="0",
        ),
        pages=[
            "Home" => "index.md",
            main_pages...,
            "Appendices" => app_pages,
            (has_refs ? ["References" => "references.md"] : [])...,
            "Code" => ["Reproducing the results" => "reproduce.md",
                [f => code_page(f) for f in CODE_FILES]...,
                "KEEP functions" => "api.md"],
        ],
    )
    unversioned(joinpath(WEB, "build"))
    return nothing
end

"""
The site is not versioned (deployed to gh-pages/jtcam/ with `versions=nothing`): drop the
`siteinfo.js` and `versions.js` scripts, which would otherwise load the version selector of
KEEP's own documentation from the gh-pages root.
"""
function unversioned(dir)
    for (root, _, files) in walkdir(dir), f in files
        endswith(f, ".html") || continue
        p = joinpath(root, f)
        s = read(p, String)
        s2 = replace(s, r"<script src=\"(\.\./)*(siteinfo|versions)\.js\"></script>" => "")
        s2 == s || write(p, s2)
    end
end

build()

if ARGS_DEPLOY
    deploydocs(;
        root=WEB, target="build", repo="$GITHUB.git", branch="gh-pages",
        dirname="jtcam", versions=nothing, push_preview=false,
        devbranch=get(ENV, "JTCAM_WEB_DEVBRANCH", "main"),
    )
end
