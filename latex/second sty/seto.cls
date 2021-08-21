%---------------------------------------------
%info
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{seto}[2021/07/10]
%% Last updated 10 July 2021
%% package name: seto.cls 
%% Copyright 2021 by Neal Yan
%
% *shhhh im stealing a friend's name lolollololol*
%
% This program can redistributed and/or modified under the terms 
% of the LaTeX Project Public License Distributed from CTAN 
% archives in directory macros/latex/base/lppl.txt; either 
% version 1 of the License, or (at your option) any later version.
%
%---------------------------------------------
\LoadClass{article}
\RequirePackage[inline]{asymptote}
\begin{asydef}
/////// Personal pens //////
pen darkgrn=RGB(13,83,76);
pen thickgrn=linewidth(1)+RGB(13,83,76);
pen redpink=RGB(229,162,155);
pen thickpink=linewidth(1)+RGB(229,162,155);
////////////////////////////
\end{asydef}
%---------------------------------------------
%options
\newif\ifnormalnums\normalnumstrue
\newif\ifheader\headertrue
\newif\ifsetogeo\setogeotrue
\newif\ifdeffonts\deffontsfalse
\DeclareOption{labelsBySect}{\normalnumsfalse}
\DeclareOption{noheader}{\headerfalse}
\DeclareOption{nogeo}{\setogeofalse}
\DeclareOption{deffonts}{\deffontstrue}
%------------------
\DeclareOption*{\ClassWarning{seto}{Unknown ‘\CurrentOption’}}
\ProcessOptions\relax
%---------------------------------------------
\ifsetogeo\RequirePackage{universal}
\else\RequirePackage[nogeo]{universal}
\fi
%colors, for mdframed
% Theme: big sur
% let's use mdframed already, its capabilities are a superset of those of tcolorbox
\RequirePackage{xstring,xcolor}
\definecolor{white}{HTML}{FFFFFF}% do not remove
\definecolor{black}{HTML}{000000}% do not remove
\definecolor{setogrn}{HTML}{0D534F}
\definecolor{probcolor}{HTML}{E5A29B} % brownish red?
\definecolor{theocolor}{HTML}{9BCAE8} % light blu
\definecolor{claimcolor}{HTML}{E6CE88} % yellow
\definecolor{remarkcolor}{HTML}{0D534F} % setogrn atm
\definecolor{quotecolor}{HTML}{91B9DC} % slightly darker blue\
%-------------------------------------------------
%fonts, text formatting
\ifdeffonts\else
\RequirePackage[cmintegrals,cmbraces]{newtxmath}
\RequirePackage{ebgaramond-maths}
\fi
\DeclareRobustCommand{\sffamily}{\lato}
\DeclareRobustCommand{\setogrn}{\color{setogrn}}
\DeclareRobustCommand{\em}{\setogrn\sffamily\bfseries}
\renewcommand{\emph}[1]{{\em #1}}
%------------------------------------------------------
% box styles
\RequirePackage[framemethod=TikZ]{mdframed}
\newenvironment{probBox}{\begin{mdframed}[
	innerleftmargin=0.35cm,innertopmargin=0.25cm,innerbottommargin=0.25cm,innerrightmargin=0.35cm,
    linecolor=probcolor,linewidth=1.5pt,backgroundcolor=probcolor!25!white,roundcorner=5pt]}
	{\end{mdframed}}
\newenvironment{claimBox}{\begin{mdframed}[
	innerleftmargin=0.35cm,innertopmargin=0.25cm,innerbottommargin=0.25cm,innerrightmargin=0.35cm,
    linecolor=claimcolor,linewidth=3pt,backgroundcolor=claimcolor!25!white,
    topline=false,bottomline=false,rightline=false,leftline=true]}
	{\end{mdframed}}
\newenvironment{thmBox}{\begin{mdframed}[
	innerleftmargin=0.35cm,innertopmargin=0.25cm,innerbottommargin=0.25cm,innerrightmargin=0.35cm,
    linecolor=theocolor,linewidth=1.5pt,backgroundcolor=theocolor!25!white,roundcorner=2pt]}
	{\end{mdframed}}
\newenvironment{remarkBox}{\begin{mdframed}[
	innerleftmargin=0.35cm,innertopmargin=0.25cm,innerbottommargin=0.25cm,innerrightmargin=0.35cm,
	linecolor=remarkcolor,linewidth=1.5pt,backgroundcolor=remarkcolor!20!white,roundcorner=2pt]}
	{\end{mdframed}}
%------------------
%cmds for convenience, because im still too lazy to define the boxes as envs
\newcounter{probnum}
\newcounter{exmpnum}
\newcounter{theonum}
\newcounter{defnum}
\newcounter{item}
\ifnormalnums
    \newcommand{\problabel}{\theprobnum}
    \newcommand{\exmplabel}{\theexmpnum}
    \newcommand{\thmlabel}{\thetheonum}
    \newcommand{\deflabel}{\thedefnum}
\else
    \newcommand{\problabel}{\thesection.\theitem}
    \newcommand{\exmplabel}{\thesection.\theitem}
    \newcommand{\thmlabel}{\thesection.\theitem}
    \newcommand{\deflabel}{\thesection.\theitem}
\fi
\newcommand{\problem}[2][]{\begin{probBox}
    {\sffamily\bfseries Problem \problabel{\mdseries\optlabel{#1}}.}\enskip#2\end{probBox}}
\newcommand{\prob}[2][]{\begin{probBox}
    {\sffamily\bfseries Problem{\mdseries\optlabel{#1}}.}\enskip#2\end{probBox}}
\newcommand{\example}[2][]{\begin{probBox}
    {\sffamily\bfseries Example \exmplabel{\mdseries\optlabel{#1}}.}\enskip#2\end{probBox}}
\newcommand{\exmp}[2][]{\begin{probBox}
    {\sffamily\bfseries Example{\mdseries\optlabel{#1}}.}\enskip#2\end{probBox}}
\newcommand{\claim}[1]{
    \begin{claimBox}{\sffamily\bfseries Claim ---\enskip}#1\end{claimBox}}
\newcommand{\lemma}[1]{
    \begin{claimBox}{\sffamily\bfseries Lemma ---\enskip}#1\end{claimBox}}
\newcommand{\theorem}[1]{\begin{thmBox}
    {\sffamily\bfseries Theorem \thmlabel.\enskip}#1\end{thmBox}}
\newcommand{\theo}[1]{\begin{thmBox}
    {\sffamily\bfseries Theorem.\enskip}#1\end{thmBox}}
\newcommand{\definition}[1]{\begin{thmBox}
    {\sffamily\bfseries Definition \deflabel.\enskip}#1\end{thmBox}}
\newcommand{\defn}[1]{\begin{thmBox}
    {\sffamily\bfseries Definition.\enskip}#1\end{thmBox}}
\newcommand{\remark}[1]{\begin{remarkBox}
    {\setogrn\sffamily\itshape Remark.\enskip}#1\end{remarkBox}}
\newcommand{\note}[1]{\begin{remarkBox}
    {\setogrn\sffamily\itshape Note.\enskip}#1\end{remarkBox}}
\ifnormalnums
    \pretocmd{\problem}{\stepcounter{probnum}}{}{}
    \pretocmd{\example}{\stepcounter{exmplabel}}{}{}
    \pretocmd{\theorem}{\stepcounter{theonum}}{}{}
    \pretocmd{\definition}{\stepcounter{defnum}}{}{}
\else
    \pretocmd{\problem}{\stepcounter{item}}{}{}
    \pretocmd{\example}{\stepcounter{item}}{}{}
    \pretocmd{\theorem}{\stepcounter{item}}{}{}
    \pretocmd{\definition}{\stepcounter{item}}{}{}
\fi
%---------------------------------------------
% headers/footers
\RequirePackage{fancyhdr}
\ifheader
    \pagestyle{fancy}
    \fancyhf{}
    \renewcommand{\headrulewidth}{0.5pt}
    \renewcommand{\footrulewidth}{0pt}
    \fancyhead[L]{\small\emph{\@title}}
    \fancyhead[R]{\small\emph{\@author, \@date}}
    \fancyfoot[C]{\thepage}
\else\fi
%---------------------------------------------
% toc
\RequirePackage[hypertexnames=false]{hyperref}
\setcounter{tocdepth}{3}
\newcommand{\toc}{
    {\em\hypersetup{
	citecolor=black,filecolor=setogrn,
	linkcolor=black,urlcolor=black}
    \tableofcontents}\newpage}
%---------------------------------------------
% title customization
\RequirePackage{titling,titlesec,setspace}
\patchcmd{\@maketitle}{\@title}{{\AlegreyaMedium\setogrn\@title}}{}{}
\patchcmd{\@maketitle}{\@author}{{\Alegreya\setogrn\@author}}{}{}
\patchcmd{\@maketitle}{\@date}{{\Alegreya\@date}}{}{}
%---------------------------------------------
% regular (sub)(sub)sections
\titleformat{\section}{\AlegreyaMedium\Large\setogrn}{\S\thesection}{1em}{}
\ifnormalnums\else
\pretocmd{\section}{\setcounter{item}{0}}{}{}
\fi
\titleformat{\subsection}{\sffamily\bfseries\large\setogrn}{\thesubsection}{1em}{}
\titleformat{\subsubsection}{\sffamily\bfseries\setogrn}{\thesubsubsection}{1em}{}
%---------------------------------------------
% misc
%-----------
% exercises
\newcounter{exenum}
\ifnormalnums
    \newcommand{\exelabel}{\theexenum}
\else
    \newcommand{\exelabel}{\thesection.\theitem}
\fi
\newcommand{\exercise}[1][]{
	\par\noindent\emph{Problem \exelabel}{\setogrn\sffamily\optlabel{#1}}\emph{.}\enskip}
\ifnormalnums
    \pretocmd{\exercise}{\stepcounter{exenum}}{}{}
\else
    \pretocmd{\exercise}{\stepcounter{item}}{}{}
\fi
%---------------------------------------------
% sol
\newenvironment{sol}[1][]{\begin{mdframed}[
    innerleftmargin=0.35cm,innertopmargin=0.15cm,innerbottommargin=0.15cm,innerrightmargin=0cm,
    topline=false,bottomline=false,rightline=false,leftline=true,
    linecolor=setogrn,linewidth=3pt,backgroundcolor=white]
    \emph{Solution}{\setogrn\sffamily\optlabel{#1}}\emph{.}\enskip}{\solqed\end{mdframed}}
%solnoteprob
\newcommand{\solnoteprob}[1]{\par\noindent{\setogrn\Alegreya\small #1}\par\noindent{\setogrn\hrule\vspace{0.2cm}}\par\noindent}
%---------------------------------------------
%quotes
\newenvironment{quoteBox}{\begin{mdframed}[
	innerleftmargin=0.35cm,innertopmargin=0.25cm,innerbottommargin=0.25cm,
	innerrightmargin=0.35cm,leftmargin=3.25in, skipabove=0cm,skipbelow=0cm,
    linewidth=2.5pt,linecolor=quotecolor,backgroundcolor=quotecolor!25!white,
    topline=false,bottomline=false,rightline=false,leftline=true]}
	{\end{mdframed}}
\renewcommand{\quote}[2]{\begin{quoteBox} {\Alegreya #1}
    {\begin{flushright}{\setogrn\Alegreya #2}\end{flushright}}\end{quoteBox}}
%---------------------------------------------
