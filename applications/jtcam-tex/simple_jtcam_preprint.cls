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
\ProvidesClass{simple_jtcam_preprint}[2020/04/16 JTCAM preprint template v1.0]

% Options for this jtcam_preprint class
\def\optbib{0}
\DeclareOption{OptBiber}{\xdef\optbib{1}}
\DeclareOption{OptSoft}{\xdef\optbib{2}}

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions*\relax
\ExecuteOptions{}
\LoadClass[11pt,a4paper]{article}


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
\RequirePackage{authblk}
\RequirePackage[english]{babel}
\RequirePackage[fleqn]{amsmath}
\RequirePackage{amsthm}
\RequirePackage{float}
\RequirePackage{algpseudocode} 
\RequirePackage{multicol}
\RequirePackage[notocbasic]{nomencl}
\RequirePackage[expansion=alltext,shrink=50]{microtype}
\RequirePackage[utf8]{inputenc}
\RequirePackage[T1]{fontenc}
\RequirePackage{siunitx}
\RequirePackage{ifmtarg}
\RequirePackage{csquotes}
\RequirePackage[ttscale=.875,oldstyle]{libertine}
\RequirePackage[libertine]{newtxmath}
\RequirePackage{geometry}
\ifnum\optbib=0
\RequirePackage[%
	backend=bibtex,
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

\RequirePackage{hyperref}
\RequirePackage[numbered,open,openlevel=1]{bookmark}
\RequirePackage[nameinlink,noabbrev]{cleveref}