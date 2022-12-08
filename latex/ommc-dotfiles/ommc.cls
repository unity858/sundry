%---------------------------------------------
%info
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{ommc}[nya!]
%% By Neal Yan, for Online Monmouth MC
% Internal use only, please!
%
% This program can redistributed and/or modified under the terms 
% of the LaTeX Project Public License Distributed from CTAN 
% archives in directory macros/latex/base/lppl.txt; either 
% version 1 of the License, or (at your option) any later version.
%
%---------------------------------------------
\newif\ifinternal\internalfalse
\DeclareOption{internal}\internaltrue
\ProcessOptions\relax

\ifinternal\LoadClass[11pt]{article}
\RequirePackage{universal}
\else\LoadClass[12pt]{exam}
\RequirePackage[exam]{universal}\fi
% universal.sty at:::
% https://raw.githubusercontent.com/tac-nayn06/sundry/main/latex/tex/universal.sty
\RequirePackage{asymptote}
\begin{asydef}
//
// Global Asymptote definitions can be put here.
//
usepackage("bm");
texpreamble("\def\V#1{\bm{#1}}");
\end{asydef}

\def\A{A}
\def\B{\V{B}}

\newcommand\overarc[1]{{%
  \setbox9=\hbox{#1}%
  \ooalign{\resizebox{\wd9}{\height}{\texttoptiebar{\phantom{A}}}\cr#1}}}

\newcommand\answer[5]{\vspace{0.2cm}$\textbf{(A)}\ #1 \qquad\textbf{(B)}\ #2 \qquad\textbf{(C)}\ #3 \qquad\textbf{(D)}\ #4 \qquad\textbf{(E)}\ #5$ \vspace{0.3cm}}

%---------------------------------------------
% PALETTE
% let's use mdframed already, its capabilities are a superset of those of tcolorbox
\RequirePackage{xstring,xcolor,lmodern} % color and font packages because lmodern > cmodern
\definecolor{main}{HTML}{ab3a42} % ommc cube maroon-- top face
\DeclareRobustCommand\sffamily\lmss
\DeclareRobustCommand\main{\color{main}}
\DeclareRobustCommand\em{\main\bfseries\sffamily}
\renewcommand\emph[1]{{\em#1}}
%------------------------------------------------------
% box styles
\RequirePackage[framemethod=TikZ]{mdframed}
\RequirePackage[most,skins]{tcolorbox}
%beamer block
\newenvironment{block}[1][Block]{\begin{tcolorbox}[enhanced,
    arc=1mm,left=.27cm,right=.3cm,top=3pt,bottom=3pt,
    colframe=main,colback=main!17!white, 
    title=#1,fonttitle={\bfseries\sffamily},boxrule=0mm,toprule=0pt,
    toptitle=1pt,bottomtitle=1pt,titlerule=0pt,drop fuzzy shadow=black]}
    {\end{tcolorbox}}

% claim-like block
\newenvironment{blocktwo}[1][BlockTwo]{\begin{mdframed}[
    innerleftmargin=0.35cm,innertopmargin=0.25cm,innerbottommargin=0.25cm,innerrightmargin=0.35cm,
    linecolor=main,linewidth=3pt,backgroundcolor=main!25!white,
    topline=false,bottomline=false,rightline=false,leftline=true]{\em#1:\enskip}}
    {\end{mdframed}}
%---------------------------------------------
% headers/footers
\ifinternal\else
\renewcommand\thepage{}
\cfoot{}\fi
%---------------------------------------------
% toc
\RequirePackage[hypertexnames=false]{hyperref}
% renaming this bc lazy
\renewcommand\a[2]{\emph{\href#1#2}}
%---------------------------------------------
% title customization
\RequirePackage{titling,titlesec,setspace}
\patchcmd\@maketitle\@title{{\bfseries\lmss\main\@title}}{}{}
\patchcmd\@maketitle\@author{{\lmss\main\@author}}{}{}
\patchcmd\@maketitle\@date{{\lmss\@date}}{}{}
%---------------------------------------------
% regular (sub)(sub)sections
\titleformat\section{\bfseries\sffamily\Large\main}
    {\includegraphics[scale=.017]{ommc-blank.png}\hspace{1mm}\thesection}{1em}{}
\titleformat\subsection{\sffamily\bfseries\large\main}
    {\includegraphics[scale=.014]{ommc-blank.png}\hspace{1mm}\thesubsection}{1em}{}
\titleformat\subsubsection{\sffamily\bfseries\main}
    {\includegraphics[scale=.011]{ommc-blank.png}\hspace{1mm}\thesubsubsection}{1em}{}
%---------------------------------------------
% for the "PROBLEMS" and "INSTRUCTIONS" embellishments
\newcommand\embellish[1]{\noindent\makebox[\linewidth]{\rule\textwidth{0.4pt}}\vspace{2mm}
    \begin{center}{\em\centering\LARGE#1}\end{center}
    \noindent\makebox[\linewidth]{\rule\textwidth{0.4pt}}}

% from seto.cls
\newcounter{exenum}
\newcommand\exelabel\theexenum
\newcommand\exercise[1][]{
    \vspace{7pt}\par\noindent\emph{Problem \exelabel}{\main\sffamily\optlabel{#1}}\emph{.}\enskip}
 \pretocmd\exercise{\stepcounter{exenum}}{}{}
%----------------------------------
