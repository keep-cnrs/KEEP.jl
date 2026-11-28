-- paper.lua: turns the pandoc AST of main.tex into Documenter pages.
--
-- Run by web/make.jl, after --citeproc, with the metadata
--   webdata: JSON file written by make.jl (labels from main.aux, code regions, TikZ figures)
--   outdir : directory where the pages (*.md) and pages.json are written
-- The pandoc output itself is discarded.
--
-- Each level-1 header starts a page. Text is rendered by pandoc's HTML writer (maths left
-- to MathJax) and wrapped in ```@raw html blocks; subsection headers are emitted as
-- Documenter headers (sidebar, anchors) and the code of figures and tables as collapsible
-- ```julia blocks. Numbers of sections, equations, figures and tables come from main.aux,
-- so that they are those of the PDF.

local data, outdir
local labels = {}   -- label -> { num = "3", kind = "equation" | "figure" | "table" | "section" | ... }
local page_of = {}  -- label or id -> page file name (e.g. "section-1.md")

-- ---------------------------------------------------------------- helpers

local function readfile(path)
  local f = assert(io.open(path, "r"))
  local s = f:read("a")
  f:close()
  return s
end

local function writefile(path, s)
  local f = assert(io.open(path, "w"))
  f:write(s)
  f:close()
end

local function html_escape(s)
  return (s:gsub("&", "&amp;"):gsub("<", "&lt;"):gsub(">", "&gt;"))
end

local function anchor_of(label) return label end

local function href(label, current)
  local p = page_of[label]
  local a = anchor_of(label)
  if p == nil then return "#" .. a end
  if p == current then return "#" .. a end
  return (p:gsub("%.md$", ".html")) .. "#" .. a
end

-- ---------------------------------------------------------------- siunitx in maths

local units = {
  meter = "m", metre = "m", second = "s", kilogram = "kg", kg = "kg", gram = "g",
  newton = "N", watt = "W", joule = "J", hertz = "Hz", radian = "rad", degree = "deg",
  percent = "%", kilo = "k", milli = "m", ampere = "A", volt = "V",
}

local function si_unit(u)
  -- \kilo\watt, \newton\meter\second\per\degree, \square\meter, \meter\squared ...
  local parts, per, pow, prefix = {}, false, nil, ""
  for cmd in u:gmatch("\\(%a+)") do
    if cmd == "per" then per = true
    elseif cmd == "square" then pow = 2
    elseif cmd == "cubic" then pow = 3
    elseif cmd == "squared" and #parts > 0 then
      local last = parts[#parts]; last.pow = (last.neg and -2 or 2)
    elseif cmd == "cubed" and #parts > 0 then
      local last = parts[#parts]; last.pow = (last.neg and -3 or 3)
    elseif cmd == "kilo" or cmd == "milli" then prefix = units[cmd]
    elseif units[cmd] then
      local p = pow or 1
      if per then p = -p end
      table.insert(parts, { sym = prefix .. units[cmd], pow = p, neg = per })
      per, pow, prefix = false, nil, ""
    end
  end
  if #parts == 0 then return u end
  local out = {}
  for _, p in ipairs(parts) do
    local s
    if p.sym == "deg" then s = "{}^\\circ"
    elseif p.sym == "%" then s = "\\%"
    else s = "\\mathrm{" .. p.sym .. "}" end
    if p.pow ~= 1 then s = s .. "^{" .. p.pow .. "}" end
    table.insert(out, s)
  end
  return table.concat(out, "\\,")
end

local function strip_braces(s) return s:sub(2, -2) end

local function siunitx(m)
  m = m:gsub("\\SI%s*(%b{})%s*(%b{})", function(n, u)
    local us = si_unit(strip_braces(u))
    if us == "{}^\\circ" or us == "\\%" then return strip_braces(n) .. us end
    return strip_braces(n) .. "\\," .. us
  end)
  m = m:gsub("\\si%s*(%b{})", function(u) return si_unit(strip_braces(u)) end)
  m = m:gsub("\\ang%s*(%b{})", function(n) return strip_braces(n) .. "^\\circ" end)
  m = m:gsub("\\num%s*(%b{})", function(n) return strip_braces(n) end)
  return m
end

-- ---------------------------------------------------------------- equations

-- Tag every numbered line with its number of the PDF (\tag{n}\label{x}), and mark the
-- others \notag: MathJax runs with tags = "none", so only explicit tags are shown.
local eqcounter = 0

local function number_line(line)
  if line:match("\\nonumber") or line:match("\\notag") then return line end
  eqcounter = eqcounter + 1
  local lab = line:match("\\label%s*{([^}]*)}")
  if lab then
    local L = labels[lab]
    local n = L and L.num or tostring(eqcounter)
    if L and tonumber(L.num) then eqcounter = tonumber(L.num) end
    return (line:gsub("\\label%s*{[^}]*}", "\\tag{" .. n .. "}"))
  end
  return line .. "\\tag{" .. eqcounter .. "}"
end

local function number_equation(m)
  local env = m:match("^%s*\\begin{(%a+%*?)}")
  if env == "equation" then
    local body = m:match("^%s*\\begin{equation}(.*)\\end{equation}%s*$")
    return "\\begin{equation}" .. number_line(body) .. "\\end{equation}"
  elseif env == "align" or env == "gather" or env == "multline" then
    local body = m:match("^%s*\\begin{" .. env .. "}(.*)\\end{" .. env .. "}%s*$")
    local lines = {}
    for line in (body .. "\\\\"):gmatch("(.-)\\\\") do table.insert(lines, line) end
    if lines[#lines]:match("^%s*$") then table.remove(lines) end
    for i, l in ipairs(lines) do lines[i] = number_line(l) end
    return "\\begin{" .. env .. "}" .. table.concat(lines, "\\\\") .. "\\end{" .. env .. "}"
  end
  return nil -- not a numbered environment
end

-- ---------------------------------------------------------------- first pass: where is each label

local function collect(blocks, file)
  local function reg(id) if id and id ~= "" then page_of[id] = file end end
  pandoc.Blocks(blocks):walk({
    Header = function(h) reg(h.identifier) end,
    Figure = function(f) reg(f.identifier) end,
    Table = function(t) reg(t.identifier) end,
    Div = function(d) reg(d.identifier) end,
    Span = function(s) reg(s.identifier) end,
    Math = function(m)
      if m.mathtype == "DisplayMath" then
        for lab in m.text:gmatch("\\label%s*{([^}]*)}") do reg(lab) end
      end
    end,
  })
end

-- ---------------------------------------------------------------- second pass: rewriting

local function prefix_caption(caption, prefix)
  local long = caption.long
  if #long > 0 and (long[1].t == "Plain" or long[1].t == "Para") then
    long[1].content:insert(1, pandoc.Space())
    long[1].content:insert(1, pandoc.Strong(pandoc.Str(prefix)))
  else
    long:insert(1, pandoc.Plain({ pandoc.Strong(pandoc.Str(prefix)) }))
  end
  caption.long = long
  return caption
end

local function rewrite(blocks, file)
  return pandoc.Blocks(blocks):walk({
    Link = function(l)
      local rtype = l.attributes["reference-type"]
      if rtype then
        local lab = l.attributes["reference"]
        local L = labels[lab]
        local n = L and L.num or ("??" .. lab)
        if L == nil then io.stderr:write("paper.lua: undefined reference " .. lab .. "\n") end
        local txt = (rtype == "eqref") and ("(" .. n .. ")") or n
        return pandoc.Link(txt, href(lab, file))
      end
      local key = l.target:match("^#(ref%-.*)$")
      if key then l.target = "references.html#" .. key; return l end
    end,
    Math = function(m)
      local t = siunitx(m.text)
      if m.mathtype == "DisplayMath" then
        local env = number_equation(t)
        if env then
          local anchors = ""
          for lab in t:gmatch("\\label%s*{([^}]*)}") do anchors = anchors .. '<a id="' .. lab .. '"></a>' end
          return pandoc.RawInline("html", anchors .. '<span class="math display">' .. html_escape(env) .. "</span>")
        end
      end
      m.text = t
      return m
    end,
    Image = function(img)
      local src = img.src:gsub("%.pdf$", ""):gsub("%.png$", "")
      img.src = src .. ".png"
      img.attributes = {}
      return img
    end,
    Figure = function(f)
      local L = labels[f.identifier]
      local sa = data.standalone[f.identifier]
      if sa then
        f.content = { pandoc.Plain({ pandoc.Image({}, sa) }) }
      end
      if L then f.caption = prefix_caption(f.caption, "Figure " .. L.num .. ":") end
      return f
    end,
    Table = function(t)
      local id = t.identifier
      return t
    end,
  })
end

-- tables: the label is on the Table, or on a Div around it
local function number_tables(blocks)
  return pandoc.Blocks(blocks):walk({
    Div = function(d)
      local L = labels[d.identifier]
      if L and L.kind == "table" then
        d.content = d.content:walk({
          Table = function(t)
            t.caption = prefix_caption(t.caption, "Table " .. L.num .. ":")
            return t
          end,
        })
        return d
      end
    end,
    Table = function(t)
      local L = labels[t.identifier]
      if L and L.kind == "table" then t.caption = prefix_caption(t.caption, "Table " .. L.num .. ":") end
      return t
    end,
  })
end

-- ---------------------------------------------------------------- rendering

local wopts = pandoc.WriterOptions({ html_math_method = "mathjax", wrap_text = "none" })

local function to_html(blocks)
  if #blocks == 0 then return "" end
  return pandoc.write(pandoc.Pandoc(blocks), "html", wopts)
end

local function code_blocks(id)
  local items = data.webcode[id]
  if not items then return "" end
  local out = {}
  table.insert(out, '```@raw html\n<details class="web-code"><summary>Code that produced this '
    .. (id:match("^tab:") and "table" or "figure") .. "</summary>\n```\n")
  for _, it in ipairs(items) do
    local where = it.file .. (it.region ~= "" and (" (region " .. it.region .. ")") or "")
    local link = it.url ~= "" and ("[" .. where .. "](" .. it.url .. ")") or ("`" .. where .. "`")
    table.insert(out, link .. ", see also [the full file](" .. it.page .. ").\n\n")
    table.insert(out, "```" .. it.lang .. "\n" .. it.code .. "\n```\n\n")
  end
  table.insert(out, "```@raw html\n</details>\n```\n")
  return table.concat(out)
end

local function block_id(b)
  if b.t == "Figure" or b.t == "Table" or b.t == "Div" then return b.identifier end
  return nil
end

local function header_text(h)
  local L = labels[h.identifier]
  local title = pandoc.utils.stringify(h.content)
  if L and L.num ~= "" then title = L.num .. " " .. title end
  return title
end

local function render_page(page)
  local out = {}
  local chunk = pandoc.Blocks({})
  local function flush()
    if #chunk > 0 then
      table.insert(out, "```@raw html\n" .. to_html(chunk) .. "\n```\n\n")
      chunk = pandoc.Blocks({})
    end
  end
  table.insert(out, "# [" .. page.title .. "](@id " .. page.id .. ")\n\n")
  for _, b in ipairs(page.blocks) do
    if b.t == "Header" then
      flush()
      local hashes = string.rep("#", math.min(b.level, 6))
      table.insert(out, hashes .. " [" .. header_text(b) .. "](@id " .. b.identifier .. ")\n\n")
    else
      chunk:insert(b)
      local id = block_id(b)
      if id and data.webcode[id] then
        flush()
        table.insert(out, code_blocks(id))
      end
    end
  end
  flush()
  return table.concat(out)
end

-- ---------------------------------------------------------------- main

function Pandoc(doc)
  data = pandoc.json.decode(readfile(pandoc.utils.stringify(doc.meta.webdata)))
  outdir = pandoc.utils.stringify(doc.meta.outdir)
  labels = data.labels

  -- split at level-1 headers; citeproc's bibliography (Div "refs") gets its own page
  local pages, cur = {}, nil
  local front = pandoc.Blocks({})
  local refs = nil
  for _, b in ipairs(doc.blocks) do
    if b.t == "Div" and b.identifier == "refs" then
      refs = b
    elseif b.t == "Header" and b.identifier == "bibliography" then
      -- citeproc's own section title: the references get a page of their own
    elseif b.t == "Header" and b.level == 1 then
      local L = labels[b.identifier]
      local num = L and L.num or ""
      local file
      if num == "" then file = (b.identifier:gsub("^sec:", "")) .. ".md"
      elseif num:match("^%d+$") then file = "section-" .. num .. ".md"
      else file = "appendix-" .. num .. ".md" end
      local title = pandoc.utils.stringify(b.content)
      if num ~= "" then title = num .. " " .. title end
      cur = { id = b.identifier, file = file, title = title, blocks = pandoc.Blocks({}),
              appendix = (num ~= "" and not num:match("^%d+$")) }
      table.insert(pages, cur)
    elseif cur then
      cur.blocks:insert(b)
    else
      front:insert(b)
    end
  end
  if refs then
    table.insert(pages, { id = "references", file = "references.md", title = "References",
                          blocks = pandoc.Blocks({ refs }) })
  end

  for _, p in ipairs(pages) do
    page_of[p.id] = p.file
    collect(p.blocks, p.file)
  end
  -- citeproc entries: ids "ref-<key>" live on the references page
  for _, p in ipairs(pages) do
    p.blocks = number_tables(rewrite(p.blocks, p.file))
  end

  local list = {}
  for _, p in ipairs(pages) do
    writefile(outdir .. "/" .. p.file, render_page(p))
    table.insert(list, { title = p.title, file = p.file, appendix = p.appendix or false })
  end
  writefile(outdir .. "/pages.json", pandoc.json.encode(list))
  writefile(outdir .. "/front.html", to_html(front))
  return doc
end
