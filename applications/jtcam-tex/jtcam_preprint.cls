%     This program is free software: you can redistribute it and/or modify
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.

%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.

%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{jtcam_preprint}[2020/04/16 JTCAM preprint template v1.0]

% Options for this jtcam_preprint class
\def\optbib{0}
\DeclareOption{OptBiber}{\xdef\optbib{1}}
\DeclareOption{OptSoft}{\xdef\optbib{2}}

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions*\relax
\ExecuteOptions{}
\LoadClass[11pt,a4paper]{article}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some of the needed packages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\RequirePackage{adjustbox}
\RequirePackage{etoolbox,xpatch}
\RequirePackage[table]{xcolor}
\RequirePackage{graphicx}
\RequirePackage{titlesec}
\RequirePackage{fancyhdr}
\RequirePackage[style]{abstract}
\RequirePackage[marginal]{footmisc}
\RequirePackage{subcaption}
\RequirePackage{enumitem}
\setlist{nosep}
\RequirePackage{authblk}
\RequirePackage[english]{babel}
\RequirePackage[fleqn]{amsmath}
\RequirePackage{amsthm}
\RequirePackage{float}
\newfloat{Algorithm}{ht}{lop}
\RequirePackage{algpseudocode} 
\RequirePackage{multicol}
\RequirePackage[notocbasic]{nomencl}
\RequirePackage[expansion=alltext,shrink=50]{microtype}
\RequirePackage[utf8]{inputenc}
\RequirePackage[T1]{fontenc}
\RequirePackage{siunitx}
\RequirePackage{ifmtarg}
\RequirePackage{csquotes}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nomenclature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%\xpatchcmd{\thenomenclature}
%  {\section*{\nomname}}
%  {}
%  {}{}
  
\xpatchcmd{\thenomenclature}{%
  \section*{\nomname}% Look for `\section*... etc.
}{% Replace it by 'nothing'
}{\typeout{Success}}{\typeout{Failure}}

\makenomenclature

\renewcommand{\nomgroup}[1]%
{%
    \ifthenelse{\equal{#1}{A}}%
    {\item[]\hspace*{-\leftmargin}%
    {\sffamily\textbf{Roman symbols}}\nolinebreak}%
    {%
        \ifthenelse{\equal{#1}{B}}%
        {\vspace{2\parsep}\item[]\hspace*{-\leftmargin}%
        {\sffamily\textbf{Greek symbols}}\nolinebreak}%
        {%
            \ifthenelse{\equal{#1}{C}}%
			{\vspace{2\parsep}\item[]\hspace*{-\leftmargin}%
			{\sffamily\textbf{Abbreviations}}\nolinebreak}%
			{%
				\ifthenelse{\equal{#1}{D}}%
				 {\vspace{2\parsep}\item[]\hspace*{-\leftmargin}%
				{\sffamily\textbf{Maths}}\nolinebreak}
				{}%
			}%
        }%
    }%
}%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Journal environments and definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\def\@keywordsname{Keywords}
\newcommand{\keywords}[1]{\def\@keywords{#1}}
\newcommand{\runningtitle}[1]{\def\@runningtitle{#1}}  
\newcommand{\runningauthor}[1]{\def\@runningauthor{#1}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fonts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\RequirePackage[ttscale=.875,oldstyle]{libertine}
\renewcommand*\libertine@figurestyle{LF}
\RequirePackage[libertine]{newtxmath}
\renewcommand*\libertine@figurestyle{OsF}
\input glyphtounicode
\pdfgentounicode=1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Page set-up/margins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\RequirePackage{geometry}%
\geometry{%
	left=32mm,
	top=17.5mm,
	bottom=17.5mm,
	right=32mm,
	footskip=10.5mm,
	footnotesep=\baselineskip,
	headheight=1cm,
	headsep=5.5mm,
	a4paper,
	}% Jan 10, 2022	
	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maths content and theorem environments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
			
% Math indent
\setlength{\mathindent}{\parindent}

% After list indents: https://tex.stackexchange.com/q/5526
\def\@endtheorem{\endtrivlist}%

% Alignment problem in align env
\patchcmd{\start@aligned}{\null\,}{\null}{}{}%

\newtheoremstyle{exampstyle}
  {1.25\topsep} 	% Space above
  {1.25\topsep} 	% Space below
  {} 			% Body font
  {} 			% Indent amount
  {\bfseries} 	% Theorem head font
  {} 			% Punctuation after theorem head
  {.5em} 		% Space after theorem head
  {\bfseries\thmname{#1}\thmnumber{ #2}\@ifmtarg{#3}{\relax}{\space---\space\thmnote{#3}}}			% Theorem head spec (can be left empty, meaning `normal')

\newtheoremstyle{proofstyle}
  {1.25\topsep} 	% Space above
  {1.25\topsep} 	% Space below
  {} 			% Body font
  {} 			% Indent amount
  {\bfseries} 	% Theorem head font
  {} 			% Punctuation after theorem head
  {.5em} 		% Space after theorem head
  {\bfseries\thmname{#1}\thmnumber{ #2}\@ifmtarg{#3}{\relax}{\space---\space\thmnote{#3}}}			% Theorem head spec (can be left empty, meaning `normal')

\theoremstyle{exampstyle}
\newtheorem{theorem}{Theorem}
\newtheorem{corollary}{Corollary}
\newtheorem{lemma}{Lemma}
\newtheorem*{Sremark}{Remark}
\newtheorem{remark}{Remark}
\newtheorem{example}{Example}
\newtheorem{proposition}{Proposition}
\newtheorem{definition}{Definition}
\theoremstyle{proofstyle}
\newtheorem*{myproof}{Proof}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List of authors and affiliations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\renewcommand{\Authfont}{\bfseries}
\renewcommand{\Affilfont}{\footnotesize\normalfont}
\renewcommand\AB@authnote[1]{\textsuperscript{#1}}
\renewcommand\AB@affilnote[1]{\textsuperscript{#1}\,}
\setlength{\affilsep}{.75\baselineskip}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORCID Icon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newbox{\myorcidaffilbox}
\sbox{\myorcidaffilbox}{\includegraphics[height=1.7ex]{ORCIDlogo}}
\newcommand{\orcid}[1]{\kern1pt\href{https://orcid.org/#1}{\usebox{\myorcidaffilbox}}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title page
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\renewcommand{\@maketitle}{%%
  \mbox{}%\vskip-2.25cm%
  \noindent%\hskip-5.25cm%
  {\setlength\parskip{\z@}\noindent}
  \hbox{}\vskip\dimexpr.3cm%
  \raggedright%
  \bookmark[named=FirstPage]{\@title}%
  {\LARGE\raggedright\bfseries\sffamily\@title\par}%
  \vskip1.5\baselineskip%
  \@author\par%
  \vskip-\baselineskip%
  \normalfont\normalsize%
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Abstract
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\renewcommand{\abstitlestyle}[1]{\relax}
\setlength{\absparsep}{\z@}
\setlength{\abstitleskip}{25pt}
\setlength{\absleftindent}{\z@}
\setlength{\absrightindent}{\absleftindent}

\appto{\abstract}{%
\pdfbookmark[section]{Abstract}{abstract}%
  \setlength{\parindent}{\z@}%
  \noindent\ignorespaces%
}
\appto{\endabstract}{%
  \noindent\footnotesize\@keywordsname\quad\@keywords\\[10pt]%
  {\rule{\textwidth}{1\p@}}\par%
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sectioning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\setcounter{secnumdepth}{3}
\titlespacing*{\section}{\z@}{\baselineskip}{.3em}
\titlespacing*{\subsection}{\z@}{.7\baselineskip}{.1em}
\titlespacing*{\subsubsection}{\z@}{.5\baselineskip}{\z@}
\titlespacing*{\paragraph}{\z@}{.4\baselineskip}{1em}
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Colored tables and algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\renewcommand{\arraystretch}{1.25}
\BeforeBeginEnvironment{tabular}{%
	\centering\small\begin{adjustbox}{bgcolor=black!6}%
	}%
\AfterEndEnvironment{tabular}{%
	\end{adjustbox}%
	}%
\BeforeBeginEnvironment{tabular*}{%
	\centering\small\begin{adjustbox}{bgcolor=black!6}%
	}%
\AfterEndEnvironment{tabular*}{%
	\end{adjustbox}%
	}%
\BeforeBeginEnvironment{algorithmic}{%
	\centering\small%
	\begin{adjustbox}{cfbox=black!6 .5pt 5pt,bgcolor=black!6}%%
	\begin{minipage}[t]{\dimexpr\linewidth-11pt}%
	}%
\AfterEndEnvironment{algorithmic}{%
	\end{minipage}%
	\end{adjustbox}%
	}%	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Page headers and footers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\fancyhf{}
\renewcommand{\footrulewidth}{.5\p@}
\renewcommand{\headrulewidth}{.5\p@}

\fancyhead[L]{%
	\footnotesize\@runningauthor\hfill\@runningtitle%
	}%
	
% all other pages
\def\version{}
\pagestyle{plain}%
\fancyfoot[L]{%
	\footnotesize\version\hfill\thepage%
}%

% first page
\pagestyle{fancy}%
\fancypagestyle{plain}{%
	\fancyhf{}
	\renewcommand{\headrule}{}
	\renewcommand{\footrulewidth}{.5\p@}
	\fancyfoot[L]{%
		\footnotesize\version\hfill\thepage}%
	}%
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Captions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\captionsetup*{%
	format=plain,%
	singlelinecheck=false,%
	font={small},%
	labelfont={small,normalfont,bf,sf},%
	%labelsep=none,%
	}
	
\captionsetup*[sub]{%
    justification=centering,
	labelfont={footnotesize,bf,sf},
	%format=plain,%
	labelformat=simple
	}
	
\renewcommand\thesubfigure{(\alph{subfigure})}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Bibliography
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\ifnum\optbib=0
\RequirePackage[%
	backend=bibtex,
	firstinits=true,
	isbn=false,%
	url=false,%
	dashed=false,			  
	style=authoryear,
	maxcitenames=1,
	maxbibnames=99,
	]{biblatex}
\fi

\ifnum\optbib=1
\RequirePackage[%
	backend=biber,
	firstinits=true,
	isbn=false,%
	url=false,%
	dashed=false,			  
	style=authoryear,
	maxcitenames=1,
	maxbibnames=99,
	]{biblatex}
\fi

\ifnum\optbib=2
\RequirePackage[%
	backend=biber,
	firstinits=true,
	isbn=false,%
    datamodel=software,
	dashed=false,			  
	style=authoryear,
	maxcitenames=1,
	maxbibnames=99,
	]{biblatex}
\usepackage{software-biblatex}

% software-biblatex specific bibliography options
\ExecuteBibliographyOptions{
	halid=true,
	swhid=true,
	swlabels=true,
	vcs=true,
	license=false}
\fi

% oai/HAL eprint type
\DeclareFieldFormat{eprint:oai}{%
\textsc{hal}:~\href{https://hal.archives-ouvertes.fr/#1}%
	{\(\langle\)\nolinkurl{#1}\(\rangle\)}
}

% other biblatex options
\DeclareFieldFormat{postnote}{\printtext[bibhyperref]{#1}}	

\renewbibmacro{in:}{}
\renewcommand{\finentrypunct}{}

\DeclareCiteCommand{\cite}
  {\usebibmacro{prenote}}
  {\usebibmacro{citeindex}%
   \printtext[bibhyperref]{\usebibmacro{cite}}}
  {\multicitedelim}
  {\usebibmacro{postnote}}

\DeclareCiteCommand*{\cite}
  {\usebibmacro{prenote}}
  {\usebibmacro{citeindex}%
   \printtext[bibhyperref]{\usebibmacro{citeyear}}}
  {\multicitedelim}
  {\usebibmacro{postnote}}

\DeclareCiteCommand{\parencite}[\mkbibparens]
  {\usebibmacro{prenote}}
  {\usebibmacro{citeindex}%
    \printtext[bibhyperref]{\usebibmacro{cite}}}
  {\multicitedelim}
  {\usebibmacro{postnote}}

\DeclareCiteCommand*{\parencite}[\mkbibparens]
  {\usebibmacro{prenote}}
  {\usebibmacro{citeindex}%
    \printtext[bibhyperref]{\usebibmacro{citeyear}}}
  {\multicitedelim}
  {\printtext[bibhyperref]{\usebibmacro{postnote}}}

\DeclareCiteCommand{\footcite}[\mkbibfootnote]
  {\usebibmacro{prenote}}
  {\usebibmacro{citeindex}%
  \printtext[bibhyperref]{ \usebibmacro{cite}}}
  {\multicitedelim}
  {\usebibmacro{postnote}}

\DeclareCiteCommand{\footcitetext}[\mkbibfootnotetext]
  {\usebibmacro{prenote}}
  {\usebibmacro{citeindex}%
   \printtext[bibhyperref]{\usebibmacro{cite}}}
  {\multicitedelim}
  {\usebibmacro{postnote}}

\DeclareCiteCommand{\textcite}
  {\boolfalse{cbx:parens}}
  {\usebibmacro{citeindex}%
   \printtext[bibhyperref]{\usebibmacro{textcite}}}
  {\ifbool{cbx:parens}
     {\bibcloseparen\global\boolfalse{cbx:parens}}
     {}%
   \multicitedelim}
  {\usebibmacro{textcite:postnote}}
  

% https://tex.stackexchange.com/a/25972/371
\DeclareCiteCommand{\textcite}
 {\boolfalse{cbx:parens}}
  {\usebibmacro{citeindex}%
    \printtext[bibhyperref]{\printnames{labelname}%
      \printtext{ (\printfield{year}\printtext{)}}}}
 {\ifbool{cbx:parens}
  {\bibcloseparen\global\boolfalse{cbx:parens}}
  {}%
 \multicitedelim}
{\usebibmacro{textcite:postnote}}

\DeclareFieldFormat{biblabeldate}{#1}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hyperlinks and bookmarks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\RequirePackage{hyperref}
\hypersetup{colorlinks=true,allcolors=blue!50!black}
\RequirePackage[numbered,open,openlevel=1]{bookmark}
\urlstyle{same}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Labels and referencing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\RequirePackage[nameinlink,noabbrev]{cleveref}
\crefname{Algorithm}{Algorithm}{Algorithms}
\crefname{proposition}{Proposition}{Propositions}
\Crefname{proposition}{Proposition}{Propositions}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% At BeginDocument
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\AtBeginDocument{
\flushbottom\maketitle	% automatically add maketitle contents
}
