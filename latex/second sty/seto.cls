%---------------------------------------------
%info
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{seto}[2021/09/06]
%% Last updated 21 Nov '22
%% package name: seto.cls 
%% Copyright 2021 by Neal Yan
%
% This program can redistributed and/or modified under the terms 
% of the LaTeX Project Public License Distributed from CTAN 
% archives in directory macros/latex/base/lppl.txt; either 
% version 1 of the License, or (at your option) any later version.
%
%---------------------------------------------
\LoadClass[11pt]{article}
\RequirePackage{asymptote}
%---------------------------------------------
%options
\newif\ifnormalnums\normalnumstrue
\newif\ifheader\headertrue
\newif\ifsetogeo\setogeotrue
\newif\ifdeffonts\deffontsfalse
\DeclareOption{labelsBySect}\normalnumsfalse
\DeclareOption{noheader}\headerfalse
\DeclareOption{nogeo}\setogeofalse
\DeclareOption{deffonts}\deffontstrue
%------------------
\DeclareOption*{\PassOptionsToClass\CurrentOption{article}}
\ProcessOptions\relax
%---------------------------------------------
\ifsetogeo\RequirePackage{universal}
\else\RequirePackage[nogeo]{universal}
\fi
% PALETTE
% let's use mdframed already, its capabilities are a superset of those of tcolorbox
\RequirePackage{xstring,xcolor}
\definecolor{white}{HTML}{FFFFFF}% do not remove
\definecolor{black}{HTML}{000000}% do not remove
\definecolor{maincolor}{HTML}{0D534F}
\definecolor{exmpcolor}{HTML}{0d534f} % brownish red?
\definecolor{theocolor}{HTML}{9BCAE8} % light blu
\definecolor{claimcolor}{HTML}{E6CE88} % yellow
\definecolor{remarkcolor}{HTML}{0D534F} % maincolor atm
\definecolor{quotecolor}{HTML}{91B9DC} % slightly darker blue
\definecolor{trunk}{HTML}{4f282d}
%-------------------------------------------------
%fonts, text formatting
\ifdeffonts\else
\RequirePackage[T1]{fontenc}
% \RequirePackage{baskervillef}
%\RequirePackage[varqu,varl,var0]{inconsolata}
%\RequirePackage[scale=.95,type1]{cabin}
%\RequirePackage[baskerville,vvarbb]{newtxmath}
%\RequirePackage[cal=boondoxo]{mathalfa}
\usepackage[cmintegrals,cmbraces]{newtxmath}
\RequirePackage{ebgaramond-maths}
\fi
\DeclareRobustCommand\sffamily\lato
\DeclareRobustCommand\maincolor{\color{maincolor}}
\DeclareRobustCommand\em{\maincolor\sffamily\bfseries}
\renewcommand\emph[1]{{\em #1}}
%------------------------------------------------------
% box styles
\RequirePackage[framemethod=TikZ]{mdframed}
\RequirePackage[most,skins]{tcolorbox}
\newenvironment{exmpBox}{\begin{mdframed}[
    innerleftmargin=0.35cm,innertopmargin=0.25cm,innerbottommargin=0.25cm,innerrightmargin=0.35cm,
    linecolor=exmpcolor,linewidth=1.5pt,backgroundcolor=exmpcolor!25!white,roundcorner=5pt]}
    {\end{mdframed}}
\newenvironment{claimBox}{\begin{mdframed}[
    innerleftmargin=0.35cm,innertopmargin=0.25cm,innerbottommargin=0.25cm,innerrightmargin=0.35cm,
    linecolor=claimcolor,linewidth=3pt,backgroundcolor=claimcolor!25!white,
    topline=false,bottomline=false,rightline=false,leftline=true]}
    {\end{mdframed}}
\newenvironment{thmBox}{\begin{mdframed}[
    innerleftmargin=0.35cm,innertopmargin=0.25cm,innerbottommargin=0.25cm,innerrightmargin=0.35cm,
    linecolor=theocolor,linewidth=1.5pt,backgroundcolor=theocolor!25!white,leftline=false,rightline=false]}
    {\end{mdframed}}
\newenvironment{remarkBox}{\begin{mdframed}[
    innerleftmargin=0.35cm,innertopmargin=0.25cm,innerbottommargin=0.25cm,innerrightmargin=0.35cm,
    linecolor=remarkcolor,linewidth=1.5pt,backgroundcolor=remarkcolor!20!white,topline=false,bottomline=false]}
    {\end{mdframed}}
%------------------
%cmds for convenience, because im still too lazy to define the boxes as envs
\newcounter{item}
\newenvironment{block}[1][Block]{\begin{tcolorbox}[enhanced,
    arc=1mm,left=.27cm,right=.3cm,top=3pt,bottom=3pt,
    colframe=maincolor,colback=maincolor!17!white, 
    title=#1,fonttitle={\bfseries\sffamily},boxrule=0mm,toprule=0pt,
    toptitle=1pt,bottomtitle=1pt,titlerule=0pt,drop fuzzy shadow=black]}
    {\end{tcolorbox}}
\newenvironment{thm}[1][Theorem]{\begin{thmBox}{\sffamily\bfseries #1}\\}{\end{thmBox}}
\newenvironment{exmp}[1][Example]{\begin{exmpBox}{\sffamily\bfseries #1}\\}{\end{exmpBox}}
\newenvironment{remark}[1][Remark]{\begin{remarkBox}
	{\maincolor\sffamily\itshape #1.\enskip}}{\end{remarkBox}}
\newcommand\claim[2][Claim]
	{\begin{claimBox}{\sffamily\bfseries #1 --\enskip}#2\end{claimBox}}
%---------------------------------------------
% headers/footers
\RequirePackage{fancyhdr}
\newcommand\insertdisplaytitle\@title
\newcommand\insertdisplayauthor\@author
\newcommand\insertdisplaydate\@date
\def\displaytitle#1{\renewcommand\insertdisplaytitle{#1}}
\def\displayauthor#1{\renewcommand\insertdisplayauthor{#1}}
\def\displaydate#1{\renewcomand\insertdisplaydate{#1}}
\ifheader
    \pagestyle{fancy}
    \fancyhf{}
    \renewcommand\headrulewidth{0.5pt}
    \renewcommand\footrulewidth{0pt}
    \fancyhead[L]{{\small\emph\insertdisplaytitle}}
    \fancyhead[R]{{\small\em\insertdisplayauthor\ {\mdseries(\insertdisplaydate)}}}
    \fancyfoot[C]{\arabic{page}}
\else\fi
%---------------------------------------------
% toc
\renewcommand\thepage{{\maincolor\sffamily\arabic{page}}}
\RequirePackage[hypertexnames=false]{hyperref}
\setcounter{tocdepth}3
\newcommand\toc{{\sffamily\maincolor
    \hypersetup{
    citecolor=maincolor,filecolor=black,
    linkcolor=black,urlcolor=maincolor}
    \tableofcontents}
    \renewcommand\thepage{\arabic{page}}\newpage}
%---------------------------------------------
% title customization
\RequirePackage{titling,titlesec,setspace}
\patchcmd\@maketitle\@title{{\AlegreyaMedium\maincolor\@title}}{}{}
\patchcmd\@maketitle\@author{{\Alegreya\maincolor\@author}}{}{}
\patchcmd\@maketitle\@date{{\Alegreya\@date}}{}{}
%---------------------------------------------
% regular (sub)(sub)sections
\titleformat\section{\AlegreyaMedium\Large\maincolor}{
\begin{tikzpicture}[scale=0.035]
    \filldraw[] (0,0)--(1,0)--(1,2)--(3.5,2)--(2,4)--(3,4)--(1.5,6)--(2.5,6)--(1,8)--(2,8)--(0,11);
    \filldraw[] (0,0)--(-1,0)--(-1,2)--(-3.5,2)--(-2,4)--(-3,4)--(-1.5,6)--(-2.5,6)--(-1,8)--(-2,8)--(0,11);
\filldraw[trunk](-1,0)--(1,0)--(1,2)--(-1,2)--(-1,0);
\end{tikzpicture}
    \thesection}{1em}{}
\ifnormalnums\else
\pretocmd\section{\setcounter{item}{0}}{}{}
\fi
\titleformat\subsection{\sffamily\bfseries\large\maincolor}{\begin{tikzpicture}[scale=0.03]
    \filldraw[] (0,0)--(1,0)--(1,2)--(3.5,2)--(2,4)--(3,4)--(1.5,6)--(2.5,6)--(1,8)--(2,8)--(0,11);
    \filldraw[] (0,0)--(-1,0)--(-1,2)--(-3.5,2)--(-2,4)--(-3,4)--(-1.5,6)--(-2.5,6)--(-1,8)--(-2,8)--(0,11);
\filldraw[trunk](-1,0)--(1,0)--(1,2)--(-1,2)--(-1,0);
\end{tikzpicture}
\thesubsection}{1em}{}
\titleformat\subsubsection{\sffamily\bfseries\maincolor}{\begin{tikzpicture}[scale=0.027]
    \filldraw[] (0,0)--(1,0)--(1,2)--(3.5,2)--(2,4)--(3,4)--(1.5,6)--(2.5,6)--(1,8)--(2,8)--(0,11);
    \filldraw[] (0,0)--(-1,0)--(-1,2)--(-3.5,2)--(-2,4)--(-3,4)--(-1.5,6)--(-2.5,6)--(-1,8)--(-2,8)--(0,11);
\filldraw[trunk](-1,0)--(1,0)--(1,2)--(-1,2)--(-1,0);
\end{tikzpicture}
\thesubsubsection}{1em}{}
%---------------------------------------------
% misc
%-----------
% exercises
\newcounter{exenum}
\ifnormalnums
    \newcommand\exelabel\theexenum
\else
    \newcommand\exelabel{\thesection.\theitem}
\fi
\newcommand\exercise[1][]{
    \vspace{7pt}\par\noindent\emph{Problem \exelabel}{\maincolor\sffamily\optlabel{#1}}\emph{.}\enskip}
\ifnormalnums
    \pretocmd\exercise{\stepcounter{exenum}}{}{}
\else
    \pretocmd\exercise{\stepcounter{item}}{}{}
\fi
%---------------------------------------------
% sol
\newenvironment{sol}[1][]{\begin{mdframed}[
    innerleftmargin=0.35cm,innertopmargin=0.15cm,innerbottommargin=0.15cm,innerrightmargin=0cm,
    topline=false,bottomline=false,rightline=false,leftline=true,
    linecolor=maincolor,linewidth=3pt,backgroundcolor=white]
    \emph{Solution}{\maincolor\sffamily\optlabel{#1}}\emph{.}\enskip}{\solqed\end{mdframed}}
%solnoteprob
\newcommand\solnoteprob[1]{\par\noindent{\maincolor\Alegreya\small #1}\par\noindent{\maincolor\hrule\vspace{0.2cm}}\par\noindent}
%---------------------------------------------
%quotes
\newenvironment{quoteBox}{\begin{mdframed}[
    innerleftmargin=0.35cm,innertopmargin=0.25cm,innerbottommargin=0.25cm,innerrightmargin=0.35cm,
    leftmargin=3.25in,skipabove=0cm,skipbelow=0cm,
    linewidth=2.5pt,linecolor=quotecolor,backgroundcolor=quotecolor!25!white,
    topline=false,bottomline=false,rightline=false,leftline=true]}
    {\end{mdframed}}
\renewcommand\quote[2]{\begin{quoteBox} {\Alegreya #1}
    {\begin{flushright}{\maincolor\Alegreya #2}\end{flushright}}\end{quoteBox}}
%---------------------------------------------
