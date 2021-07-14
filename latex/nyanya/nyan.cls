%---------------------------------------------
%info
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{nyan}[2021/07/10]
%% package name: nyan.cls 
%% Copyright 2021 by Neal Yan
% 
% seriously, just use seto.cls already XD
%
% This program can redistributed and/or modified under the terms 
% of the LaTeX Project Public License Distributed from CTAN 
% archives in directory macros/latex/base/lppl.txt; either 
% version 1 of the License, or (at your option) any later version.
%
≈---------------------------------
% options
\newboolean{normalnumbering}
\newboolean{sectnumbering}
\newboolean{headers}
\newboolean{credit}
\setboolean{normalnumbering}{true}
\setboolean{sectnumbering}{false}
\setboolean{headers}{true}
\setboolean{credit}{true}
\DeclareOption{labelsBySect}{
    \setboolean{normalnumbering}{false}
    \setboolean{sectnumbering}{true}}
\DeclareOption{noheader}{\setboolean{headers}{false}}
\DeclareOption{nocredit}{\setboolean{credit}{false}}
%--------------
\DeclareOption*{\ClassWarning{nyan}{Unknown ‘\CurrentOption’}}
\ProcessOptions\relax
%---------------------------------------------
%fonts, text formatting, the default stuff
\LoadClass{article}
\RequirePackage{universal}
%---------------------------------------------
%colors
% Theme: nyan cat colors
\RequirePackage{xstring,xcolor}
\definecolor{probcolor}{HTML}{FF0090}%\problem/example color
\definecolor{white}{HTML}{FFFFFF}%dont remove this one
\definecolor{black}{HTML}{000000}% or this one
\definecolor{remarkcolor}{HTML}{FF9898}%\remark color
\definecolor{claimcolor}{HTML}{F5C596}%\claim color
\definecolor{nyanblue}{HTML}{003366}%
\definecolor{theocolor}{HTML}{FD99FD}%\example color
\DeclareRobustCommand{\nyanblue}{\color{nyanblue}}
%---------------------------------------------
% tcolorboxes
\RequirePackage[most]{tcolorbox}
\newtcolorbox{probBox}{colframe=probcolor,colback=probcolor!20!white}
\newtcolorbox{claimBox}{colframe=claimcolor,colback=claimcolor!25!white}
\newtcolorbox{theoBox}{colframe=theocolor,colback=theocolor!30!white}
\newtcolorbox{remarkBox}{colframe=remarkcolor,colback=remarkcolor!25!white}
%--------------
% commands using the tcolorboxes, im too lazy to def envs lul
\newcommand{\optlabel}[1]{\ifthenelse{\isempty{#1}}{}{ (#1)}}
\newcounter{probnum}
\newcommand{\problabel}{\ifthenelse{\boolean{normalnumbering}}{\theprobnum}{\thesection.\theprobnum}}
\newcommand{\problem}[2][]{\stepcounter{probnum}\begin{problemBox}
    {\cmss\bfseries Problem \problabel}{\cmss\optlabel(#1)}
    {\cmss\bfseries.}\enskip #2\end{problemBox}}
\newcommand{\prob}[2][]{\begin{problemBox}
    {\cmss\bfseries Problem}{\cmss\optlabel{#1}}
    {\cmss\bfseries.}\enskip #2\end{problemBox}}
\newcounter{exmpnum}
\newcommand{\exmplabel}{\ifthenelse{\boolean{normalnumbering}}{\theexmpnum}{\thesection.\theexmpnum}}
\newcommand{\example}[2][]{\stepcounter{exmpnum}\begin{problemBox}
	{\cmss\bfseries Example \exmplabel}{\cmss\optlabel{#1}}
	{\cmss\bfseries.}\enskip #2\end{problemBox}}
\newcommand{\exmp}[2][]{\begin{problemBox}
    {\cmss\bfseries Example}{\cmss\optlabel{#1}}
    {\cmss\bfseries.}\enskip #2\end{problemBox}}
\newcommand{\claim}[1]{\begin{claimBox}
	{\cmss\bfseries Claim ---\quad}#1\end{claimBox}} % evan [chen]'s style
\newcommand{\lemma}[1]{\begin{lemmaBox}
	{\cmss\bfseries Lemma ---\enskip}#1\end{lemmaBox}}
\newcommand{\remark}[1]{\begin{remarkBox}
	{\color{nyanblue}\cmss\itshape Remark.\enskip}#1\end{remarkBox}}
\newcounter{theonum}
\newcommand{\thmlabel}{\ifthenelse{\boolean{normalnumbering}}{\thetheonum}{\thesection.\thetheonum}}
\newcommand{\theorem}[1]{\stepcounter{theonum}\begin{theoBox}
	{\cmss\bfseries Theorem \thmlabel.\enskip}#1\end{theoBox}}
\newcommand{\theo}[1]{\begin{theoBox}
	{\cmss\bfseries Theorem.\enskip}#1\end{theoBox}}
\newcounter{defnum}
\newcommand{\deflabel}{\ifthenelse{\boolean{normalnumbering}}{\thedefnum}{\thesection.\thedefnum}}
\newcommand{\definition}[1]{\stepcounter{defnum}\begin{theoBox}
	{\cmss\bfseries Definition \deflabel.\enskip}#1\end{theoBox}}
\newcommand{\defn}[1]{\begin{theoBox}
	{\cmss\bfseries Definition.\enskip}#1\end{theoBox}}
\newcommand{\note}[1]{\begin{remarkBox}
    {\color{nyanblue}\cmss\itshape Note.\enskip}#1\end{remarkBox}}
\newcommand{\pro}[1]{\begin{proof}#1\end{proof}}
%---------------------------------------------
% text formatting
\DeclareRobustCommand{\sffamily}{\cmss}
\renewcommand{\emph}[1]{{\nyanblue\cmss\bfseries #1}}
%---------------------------------------------
% exercises
\newcounter{exenum}
\newcommand{\exelabel}{\ifthenelse{\boolean{sectnumbering}}{\thesection.\theexenum}{\theexenum}}
\newcommand{\exercise}[1][]{\stepcounter{exenum}
	\par\noindent\emph{Problem \exelabel}{\color{nyanblue}\cmss\optlabel{#1}}\emph{.}\enskip}
%---------------------------------------------
% headers/footers
\RequirePackage{fancyhdr}
\ifthenelse{\boolean{headers}}{
    \pagestyle{fancy}\fancyhf{}
    \renewcommand{\headrulewidth}{0.5pt}
    \renewcommand{\footrulewidth}{0pt}
    \fancyhead[L]{\emph{\@title}}
    \fancyhead[R]{\emph{\@author, \@date}}
    \fancyfoot[C]{\thepage}}{}
%---------
% toc
\setcounter{tocdepth}{3}
\newcommand{\toc}{{\hypersetup{citecolor=black,filecolor=nyanblue,linkcolor=black,urlcolor=black}
	\emph{\tableofcontents}\newpage}}
%---------------------------------------------
% headings
\RequirePackage{titling,titlesec,setspace}
% title
\ifthenelse{\boolean{credit}}{\apptocmd{\maketitle}{\begin{center}
    \textit{Produced using nyanya.cls}\end{center}}{}{}}{}% i'm too lazy to change font for \maketitle lol 
\patchcmd{\@maketitle}{\@title}{{\nyanblue\cmss\bfserues\@title}}{}{}
\patchcmd{\@maketitle}{\@author}{{\cmss\color{nyanblue}\@author}}{}{}
\patchcmd{\@maketitle}{\@date}{{\cmss\bfseries\@date}}{}{}
% sections
\titleformat{\section}{\nyanblue\Large\cmss\bfseries}{{\S\thesection}}{1em}{}
\titleformat{\subsection}{\nyanblue\cmss\large\bfseries}{{\thesubsection}}{1em}{}
\titleformat{\subsubsection}{\nyanblue\cmss\bfseries}{{\thesubsubsection}}{1em}{}
%---------------------------------------------
