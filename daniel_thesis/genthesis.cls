% Author: Daniel Celis Garza <daniel.celisgarza@materials.ox.ac.uk>
% Date of creation: 2017/03/08
% Last edit: 2017/04/25
% This class provides a generic thesis template based on the report document class.

%----------------------------- Header -----------------------------%
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{genthesis}[2017/03/08 v1.0 Generic thesis template.]
%------------------------------------------------------------------%

%---------------------- Document Definition -----------------------%
% Default document is one sided, no option needed. The document 
% option "twoside" makes the document double sided and open right 
% like a western book.

\DeclareOption{twoside}{
	\PassOptionsToClass{twoside, openright}{report}
}
%------------------------------------------------------------------%

%--------------------------- Logo ---------------------------------%
% An option must be specified because the document needs it to 
% redefine the \logo command. If an extra logo is to be added,
% make sure to the new value of \nlogo to the options of the 
% redefinition of \logo.

% Single logo.
% Control logo size by its width.
\newcommand{\nlogo}{wone}
% Control logo by its height.
\DeclareOption{honelogo}{\renewcommand{\nlogo}{hone}}

% Two logos.
% Control both logos by their width.
\DeclareOption{wtwologo}{\renewcommand{\nlogo}{wtwo}}
% Control both logos by their height.
\DeclareOption{htwologo}{\renewcommand{\nlogo}{htwo}}
%------------------------------------------------------------------%

%--------------------------- Supervisor ---------------------------%
% An option must be specified because the document needs it to 
% redefine the \supervisor command. If an extra supervisor is to 
% be added, make sure to the new value of \nspvr to the options of 
% the redefinition of \supervisor.

% One supervisor.
\newcommand{\nspvr}{one}
% Two supervisors.
\DeclareOption{twosup}{\renewcommand{\nspvr}{two}}
% Three supervisors.
\DeclareOption{threesup}{\renewcommand{\nspvr}{three}}
% Multiple supervisors.
\DeclareOption{multisup}{\renewcommand{\nspvr}{multi}}
%------------------------------------------------------------------%

%-------------------- Miscellaneous Arguments ---------------------%
% Pass all other document options to the base document class.
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{report}}
\ProcessOptions\relax
%------------------------------------------------------------------%

%----------------------- Load Report Class ------------------------%
\LoadClass[a4paper]{report}
%------------------------------------------------------------------%

%----------------------- Required Packages ------------------------%
% If then operations (\logo and \supervisor command definitions)
\RequirePackage{ifthen}
% Create and modify counters.
\RequirePackage{chngcntr}
% Include graphics.
\RequirePackage{graphicx}
%------------------------------------------------------------------%

%--------------------------- Formatting ---------------------------%
% Line spacing will not be adjusted if the natural spacing 
% cannot fit an integer number of lines in a page.
\raggedbottom
%------------------------------------------------------------------%

%------------------ Title Page Command Definition -----------------%
% University.
\newcommand{\university}[1]{\def\@university{#1}}
% University department.
\newcommand{\department}[1]{\def\@department{Department of #1}}
% College.
\newcommand{\college}[1]{\def\@college{#1 College}}
% Default submitted text.
\newcommand{\submittedtext}{A thesis submitted for the degree of}
% Degree.
\newcommand{\degree}[1]{\def\@degree{#1}}
% Date of award.
\newcommand{\degreedate}[1]{\def\@degreedate{#1}}

% Single logo by width.
\ifthenelse{\equal{\nlogo}{wone}}{
	\newcommand{\logo}[2]{\def\@logo{\includegraphics[width=#1]{#2}}}
}{}
% Double logo by width.
\ifthenelse{\equal{\nlogo}{wtwo}}{
	\newcommand{\logo}[4]{\def\@logo{\hfill\includegraphics[width=#1]{#2}\hfill\includegraphics[width=#3]{#4}\hfill}}
}{}
% Single logo by height.
\ifthenelse{\equal{\nlogo}{hone}}{
	\newcommand{\logo}[2]{\def\@logo{\includegraphics[height=#1]{#2}}}
}{}
% Double logo by height.
\ifthenelse{\equal{\nlogo}{htwo}}{
	\newcommand{\logo}[4]{\def\@logo{\hfill\includegraphics[height=#1]{#2}\hfill\includegraphics[height=#3]{#4}\hfill}}
}{}
% 1 Supervisor.
\ifthenelse{\equal{\nspvr}{one}}{
	\newcommand{\supervisor}[1]{\def\@supervisor{Supervisor: #1}}
}{}
% 2 Supervisors.
\ifthenelse{\equal{\nspvr}{two}}{
	\newcommand{\supervisor}[2]{\def\@supervisor{Supervisors: #1 \& #2}}
}{}
% 3 Supervisors.
\ifthenelse{\equal{\nspvr}{three}}{
	\newcommand{\supervisor}[3]{\def\@supervisor{Supervisors: #1, #2 \& #3}}
}{}
% Multiple Supervisors.
\ifthenelse{\equal{\nspvr}{multi}}{
	\newcommand{\supervisor}[1]{\def\@supervisor{Supervisors: #1}}
}{}
%------------------------------------------------------------------%

%---------------------- Force Single Column -----------------------%
% It is always single column even in double column documents. 
% Execute \clearpage if it's a single sided thesis, 
% \cleardoublepage if it's a double sided document.
\newenvironment{single}{%
	\@restonecolfalse
	\if@twocolumn
		\@restonecoltrue\onecolumn
	\else
		\if@openright
			\cleardoublepage
		\else
			\clearpage
		\fi
	\fi
}%
{
	\if@restonecol
		\twocolumn
	\else
		\newpage\thispagestyle{empty}
	\fi
}
%------------------------------------------------------------------%

%------------------------- New Maketitle --------------------------%
% Create title page.
\renewcommand{\maketitle}{%
	\begin{single}
		\pagenumbering{alph}
		\renewcommand{\footnotesize}{\small}
		\renewcommand{\footnoterule}{\relax}
		\thispagestyle{empty}
		\vfill
		\begin{center}
			% Title.
			{\Huge \textbf{\@title} \par}
			% Logo(s).
			{\large \vfill {\@logo \par} \vfill}
			% Author.
			{{\LARGE \@author} \par}
			{\Large \vspace*{1ex}
				% University.
				{\@university \par}
				\vspace*{1ex}
				% College.
				{\@college \par}
				\vspace*{1ex}
				% Department.
				{\@department \par}
				\vspace*{1ex}
				% Supervisor(s).
				{\@supervisor \par}
				\vfill
				% Submitted text.
				{\submittedtext \par}
				\vspace*{1ex}
				% Degree.
				{\emph{\@degree} \par}
				\vspace*{3ex}
				% Date.
				{\@degreedate}
			}
		\end{center}
		\vfill
	\end{single}
}
%------------------------------------------------------------------%

%-------------------------- Front Matter --------------------------%
\newenvironment{frontmatter}[1][\unskip]{%
	\begin{single}
		\thispagestyle{empty}
		\begin{center}
			{\LARGE \textbf{#1}}
		\end{center}
		\begin{center}
			\vfill
}
{
			\vfill
		\end{center}
	\end{single}
}
%------------------------------------------------------------------%

%------------------------- Roman Numerals -------------------------%
% Environment for roman numerals. Carries the roman number on for 
% subsequent calls. Resets the page number to start at 1 with 
% arabic numerals after it ends.
\newcounter{romancounter}
\setcounter{romancounter}{1}
\newcounter{arabiccounter}
\setcounter{arabiccounter}{1}
	
\newenvironment{romanpages}{%
	\cleardoublepage
	\setcounter{page}{\number\value{romancounter}}\renewcommand{\thepage}{\roman{page}}
}
{%
	\setcounter{arabiccounter}{\value{arabiccounter}}
	\cleardoublepage
	\setcounter{romancounter}{\number\value{page}}
	\pagenumbering{arabic}
	\renewcommand{\thepage}{\arabic{page}}\setcounter{page}{\value{arabiccounter}}
}
%------------------------------------------------------------------%