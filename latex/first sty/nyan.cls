%---------------------------------------------
%info
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{nyan}[2021/07/10]
%% class name: nyan.cls 
%% Copyright 2021 by Neal Yan
% 
% seriously, just use seto.cls already XD
%
% This program can redistributed and/or modified under the terms 
% of the LaTeX Project Public License Distributed from CTAN 
% archives in directory macros/latex/base/lppl.txt; either 
% version 1 of the License, or (at your option) any later version.
%
%---------------------------------
% options
\newif\ifnormalnums\normalnumstrue
\newif\ifheader\headertrue
\newif\ifnyangeo\nyangeotrue
\DeclareOption{labelsBySect}{\normalnumsfalse}
\DeclareOption{noheader}{\headerfalse}
\DeclareOption{nogeo}{\nyangeofalse}
%--------------
\DeclareOption*{\ClassWarning{nyan}{Unknown ‘\CurrentOption’}}
\ProcessOptions\relax
%---------------------------------------------
%fonts, text formatting, the default stuff
\LoadClass{article}
\ifnyangeo\RequirePackage{universal}\else\RequirePackage[nogeo]{universal}\fi
\DeclareRobustCommand\sffamily\cmss
\DeclareRobustCommand\em{\nyanblue\cmss\bfseries}
\renewcommand\emph[1]{{\em #1}}
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
\newtcolorbox{probBox}{arc=1mm,left=7pt,right=7pt,top=5pt,bottom=5pt,colframe=probcolor,colback=probcolor!20!white}
\newtcolorbox{claimBox}{arc=1mm,left=7pt,right=7pt,top=5pt,bottom=5pt,colframe=claimcolor,colback=claimcolor!25!white}
\newtcolorbox{thmBox}{arc=1mm,left=7pt,right=7pt,top=5pt,bottom=5pt,colframe=theocolor,colback=theocolor!30!white}
\newtcolorbox{remarkBox}{arc=1mm,left=7pt,right=7pt,top=5pt,bottom=5pt,colframe=remarkcolor,colback=remarkcolor!25!white}
%--------------
% commands using the tcolorboxes, im too lazy to def envs lul
\newcounter{probnum}
\newcounter{exmpnum}
\newcounter{theonum}
\newcounter{defnum}
\newcounter{item}
\ifnormalnums
    \def\problabel\theprobnum
    \def\exmplabel\theexmpnum
    \def\thmlabel\thetheonum
    \def\deflabel\thedefnum
\else
    \def\problabel{\thesection.\theitem}
    \def\exmplabel{\thesection.\theitem}
    \def\thmlabel{\thesection.\theitem}
    \def\deflabel{\thesection.\theitem}
\fi
\newcommand\problem[2][]{\begin{probBox}
    {\cmss\bfseries Problem \problabel{\mdseries\optlabel{#1}}.}\enskip #2\end{probBox}}
\newcommand\prob[2][]{\begin{probBox}
    {\cmss\bfseries Problem{\mdseries\optlabel{#1}}.}\enskip #2\end{probBox}}
\newcommand\example[2][]{\begin{probBox}
    {\cmss\bfseries Example \exmplabel{\mdseries\optlabel{#1}}.}\enskip #2\end{probBox}}
\newcommand\exmp[2][]{\begin{probBox}
    {\cmss\bfseries Example{\mdseries\optlabel{#1}}.}\enskip #2\end{probBox}}
\newcommand\claim[1]{\begin{claimBox}
    {\cmss\bfseries Claim ---\quad}#1\end{claimBox}} % evan [chen]'s style
\newcommand\lemma[1]{\begin{claimBox}
    {\cmss\bfseries Lemma ---\enskip}#1\end{claimBox}}
\newcommand\remark[1]{\begin{remarkBox}
    {\nyanblue\cmss\itshape Remark.\enskip}#1\end{remarkBox}}
\newcommand\theorem[1]{\begin{thmBox}
    {\cmss\bfseries Theorem \thmlabel.\enskip}#1\end{thmBox}}
\newcommand\theo[1]{\begin{thmBox}
    {\cmss\bfseries Theorem.\enskip}#1\end{thmBox}}
\newcommand\definition[1]{\begin{thmBox}
    {\cmss\bfseries Definition \deflabel.\enskip}#1\end{thmBox}}
\newcommand\defn[1]{\begin{thmBox}
    {\cmss\bfseries Definition.\enskip}#1\end{thmBox}}
\newcommand\note[1]{\begin{remarkBox}
    {\nyanblue\cmss\itshape Note.\enskip}#1\end{remarkBox}}
\ifnormalnums
    \pretocmd\problem{\stepcounter{probnum}}{}{}
    \pretocmd\example{\stepcounter{exmplabel}}{}{}
    \pretocmd\theorem{\stepcounter{theonum}}{}{}
    \pretocmd\definition{\stepcounter{defnum}}{}{}
\else
    \pretocmd\problem{\stepcounter{item}}{}{}
    \pretocmd\example{\stepcounter{item}}{}{}
    \pretocmd\theorem{\stepcounter{item}}{}{}
    \pretocmd\definition{\stepcounter{item}}{}{}
\fi
%---------------------------------------------
% exercises
\newcounter{exenum}
\ifnormalnums
    \newcommand\exelabel\theexenum
\else
    \newcommand\exelabel{\thesection.\theitem}
\fi
\newcommand\exercise[1][]{
	\par\noindent\emph{Problem \exelabel}{\nyanblue\cmss\optlabel{#1}}\emph{.}\enskip}
\ifnormalnums
    \pretocmd\exercise{\stepcounter{exenum}}{}{}
\else
    \pretocmd\exercise{\stepcounter{item}}{}{}
\fi
%---------------------------------------------
% headers/footers
\RequirePackage{fancyhdr}
\ifheader
    \pagestyle{fancy}
    \fancyhf{}
    \renewcommand\headrulewidth{0.5pt}
    \renewcommand\footrulewidth{0pt}
    \fancyhead[L]{\small\emph{\@title}}
    \fancyhead[R]{\small\emph{\@author, \@date}}
    \fancyfoot[C]\thepage
\else\fi
%---------
% toc
\RequirePackage[hypertexnames=false]{hyperref}
\setcounter{tocdepth}{3}
\newcommand\toc{
	{\hypersetup{
	citecolor=black,filecolor=nyanblue,
	linkcolor=black,urlcolor=black}
\emph{\tableofcontents}\newpage}}
%---------------------------------------------
% headings
\RequirePackage{titling,titlesec,setspace}
% title
\patchcmd{\@maketitle}{\@title}{{\nyanblue\cmss\bfseries\@title}}{}{}
\patchcmd{\@maketitle}{\@author}{{\nyanblue\cmss\@author}}{}{}
\patchcmd{\@maketitle}{\@date}{{\cmss\@date}}{}{}
% sections
\titleformat\section{\nyanblue\Large\cmss\bfseries}{{\S\thesection}}{1em}{}
\ifnormalnums\else
\pretocmd\section{\setcounter{item}{0}}{}{}
\fi
\titleformat\subsection{\nyanblue\cmss\large\bfseries}{{\thesubsection}}{1em}{}
\titleformat\subsubsection{\nyanblue\cmss\bfseries}{{\thesubsubsection}}{1em}{}
%---------------------------------------------
% sol
\newenvironment{sol}[1][]{\vspace{0.1cm}\\\emph{Solution}{\nyanblue\sffamily\optlabel{#1}}\emph{.}\enskip}{\solqed}
%---------------------------------------------
