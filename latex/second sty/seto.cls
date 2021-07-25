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
\RequirePackage{universal} % my own
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
\RequirePackage{xifthen}
\newboolean{normalnumbering}
\newboolean{sectnumbering}
\newboolean{headers}
\setboolean{normalnumbering}{true}
\setboolean{sectnumbering}{false}
\setboolean{headers}{true}
\DeclareOption{labelsBySect}{
    \setboolean{normalnumbering}{false}
    \setboolean{sectnumbering}{true}}
\DeclareOption{noheader}{\setboolean{headers}{false}}
%------------------
\DeclareOption*{\ClassWarning{seto}{Unknown ‘\CurrentOption’}}
\ProcessOptions\relax
%---------------------------------------------
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
\RequirePackage[cmintegrals,cmbraces]{newtxmath}
\RequirePackage{ebgaramond-maths}
\DeclareRobustCommand{\sffamily}{\lato}
\DeclareRobustCommand{\setogrn}{\color{setogrn}}
\renewcommand{\emph}[1]{{\setogrn\sffamily\bfseries #1}}
%------------------------------------------------------
% box styles
\RequirePackage[framemethod=TikZ]{mdframed}
\newenvironment{probBox}{\begin{mdframed}[
	innerleftmargin=0.4cm,innertopmargin=0.25cm,innerbottommargin=0.25cm,innerrightmargin=0.4cm,
    linecolor=probcolor,linewidth=1.5pt,backgroundcolor=probcolor!25!white,roundcorner=5pt]}
	{\end{mdframed}}
\newenvironment{claimBox}{\begin{mdframed}[
	innerleftmargin=0.4cm,innertopmargin=0.25cm,innerbottommargin=0.25cm,innerrightmargin=0.4cm,
    linecolor=claimcolor,linewidth=3pt,backgroundcolor=claimcolor!25!white,
    topline=false,bottomline=false,rightline=false,leftline=true]}
	{\end{mdframed}}
\newenvironment{thmBox}{\begin{mdframed}[
	innerleftmargin=0.4cm,innertopmargin=0.25cm,innerbottommargin=0.25cm,innerrightmargin=0.4cm,
    linecolor=theocolor,linewidth=1.5pt,backgroundcolor=theocolor!25!white,roundcorner=2pt]}
	{\end{mdframed}}
\newenvironment{remarkBox}{\begin{mdframed}[
	innerleftmargin=0.4cm,innertopmargin=0.25cm,innerbottommargin=0.25cm,innerrightmargin=0.4cm,
	linecolor=remarkcolor,linewidth=1.5pt,backgroundcolor=remarkcolor!20!white,roundcorner=2pt]}
	{\end{mdframed}}
%------------------
%cmds for convenience, because im still too lazy to define the boxes as envs
\newcommand{\optlabel}[1]{\ifthenelse{\isempty{#1}}{}{ (#1)}}
\newcounter{probnum}
\newcommand{\problabel}{\ifthenelse{\boolean{normalnumbering}}{\theprobnum}{\thesection.\theprobnum}}
\newcommand{\problem}[2][]{\stepcounter{probnum}\begin{probBox}
    {\sffamily\bfseries Problem \problabel}{\sffamily\optlabel{#1}}{\sffamily\bfseries.}\enskip#2\end{probBox}}
\newcommand{\prob}[2][]{\begin{probBox}
    {\sffamily\bfseries Problem}{\sffamily\optlabel{#1}}{\sffamily\bfseries.}\enskip#2\end{probBox}}
\newcounter{exmpnum}
\newcommand{\exmplabel}{\ifthenelse{\boolean{normalnumbering}}{\theexmpnum}{\thesection.\theexmpnum}}
\newcommand{\example}[2][]{\stepcounter{exmpnum}\begin{probBox}
    {\sffamily\bfseries Example \exmplabel}{\sffamily\optlabel{#1}}{\sffamily\bfseries.}\enskip#2\end{probBox}}
\newcommand{\exmp}[2][]{\begin{probBox}
    {\sffamily\bfseries Example}{\sffamily\optlabel{#1}}{\sffamily\bfseries.}\enskip#2\end{probBox}}
\newcommand{\claim}[1]{
    \begin{claimBox}{\sffamily\bfseries Claim ---\enskip}#1\end{claimBox}}
\newcommand{\lemma}[1]{
    \begin{claimBox}{\sffamily\bfseries Lemma ---\enskip}#1\end{claimBox}}
\newcounter{theonum}
\newcommand{\thmlabel}{\ifthenelse{\boolean{normalnumbering}}{\thetheonum}{\thesection.\thetheonum}}
\newcommand{\theorem}[1]{\stepcounter{theonum}\begin{thmBox}
    {\sffamily\bfseries Theorem \thmlabel.\enskip}#1\end{thmBox}}
\newcommand{\theo}[1]{\begin{thmBox}
    {\sffamily\bfseries Theorem.\enskip}#1\end{thmBox}}
\newcounter{defnum}
\newcommand{\deflabel}{\ifthenelse{\boolean{normalnumbering}}{\thedefnum}{\thesection.\thedefnum}}
\newcommand{\definition}[1]{\stepcounter{defnum}\begin{thmBox}
    {\sffamily\bfseries Definition \deflabel.\enskip}#1\end{thmBox}}
\newcommand{\defn}[1]{\begin{thmBox}
    {\sffamily\bfseries Definition.\enskip}#1\end{thmBox}}
\newcommand{\remark}[1]{\begin{remarkBox}
    {\setogrn\sffamily\itshape Remark.\enskip}#1
    \end{remarkBox}}
\newcommand{\note}[1]{\begin{remarkBox}
    {\setogrn\sffamily\itshape Note.\enskip}#1
    \end{remarkBox}}
%---------------------------------------------
% headers/footers
\RequirePackage{fancyhdr}
\ifthenelse{\boolean{headers}}{
    \pagestyle{fancy}
    \fancyhf{}
    \renewcommand{\headrulewidth}{0.5pt}
    \renewcommand{\footrulewidth}{0pt}
    \fancyhead[L]{\small\emph{\@title}}
    \fancyhead[R]{\small\emph{\@author, \@date}}
    \fancyfoot[C]{\thepage}}{}
%---------
% toc
\RequirePackage[hypertexnames=false]{hyperref}
\setcounter{tocdepth}{3}
\newcommand{\toc}{
    {\hypersetup{
	citecolor=black,filecolor=setogrn,
	linkcolor=black,urlcolor=black}
    \emph{\tableofcontents}\newpage}}
%---------------------------------------------
% title customization
\RequirePackage{titling,titlesec,setspace}
\patchcmd{\@maketitle}{\@title}{{\AlegreyaMedium\setogrn\@title}}{}{}
\patchcmd{\@maketitle}{\@author}{{\Alegreya\setogrn\@author}}{}{}
\patchcmd{\@maketitle}{\@date}{{\Alegreya\@date}}{}{}
%---------------------------------------------
% regular (sub)(sub)sections
\titleformat{\section}{\AlegreyaMedium\Large\setogrn}{\S\thesection}{1em}{}
\apptocmd{\@section}{\ifthenelse{\boolean{sectnumbering}}{
    \setcounter{probnum}{0}\setcounter{exmpnum}{0}
    \setcounter{theonum}{0}\setcounter{defnum}{0}
    \setcounter{equation}{0}}{}}{}{}
\titleformat{\subsection}{\sffamily\bfseries\large\setogrn}{\thesubsection}{1em}{}
\titleformat{\subsubsection}{\sffamily\bfseries\setogrn}{\thesubsubsection}{1em}{}
%---------------------------------------------
% misc
%-----------
% exercises
\newcounter{exenum}
\newcommand{\exelabel}{\ifthenelse{\boolean{sectnumbering}}{\thesection.\theexenum}{\theexenum}}
\newcommand{\exercise}[1][]{\stepcounter{exenum}
	\par\noindent\emph{Problem \exelabel}{\setogrn\sffamily\optlabel{#1}}\emph{.}\enskip}
%---------------------------------------------
% sol
\newenvironment{sol}[1][]{\begin{mdframed}[
    innerleftmargin=0.4cm,innertopmargin=0.15cm,innerbottommargin=0.15cm,innerrightmargin=0cm,
    topline=false,bottomline=false,rightline=false,leftline=true,
    linecolor=setogrn,linewidth=3pt,backgroundcolor=white]
    \emph{Solution}{\setogrn\sffamily\optlabel{#1}}\emph{.}\enskip}{\solqed\end{mdframed}}
%---------------------------------------------
%quotes
\newenvironment{quoteBox}{\begin{mdframed}[
	innerleftmargin=0.4cm,innertopmargin=0.25cm,innerbottommargin=0.25cm,
	innerrightmargin=0.4cm,leftmargin=3.25in, skipabove=0cm,skipbelow=0cm,
    linewidth=2.5pt,linecolor=quotecolor,backgroundcolor=quotecolor!25!white,
    topline=false,bottomline=false,rightline=false,leftline=true]}
	{\end{mdframed}}
\renewcommand{\quote}[2]{\begin{quoteBox} {\Alegreya #1}
    {\begin{flushright}{\setogrn\Alegreya #2}\end{flushright}}
    \end{quoteBox}}
%---------------------------------------------
