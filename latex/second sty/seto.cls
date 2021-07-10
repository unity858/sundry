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
\newboolean{normalnumbering}
\newboolean{sectnumbering}
\newboolean{headers}
\newboolean{credit}
\newboolean{fancyfont}
\setboolean{normalnumbering}{true}
\setboolean{sectnumbering}{false}
\setboolean{headers}{true}
\setboolean{credit}{true}
\setboolean{fancyfont}{true}
\DeclareOption{labelsBySect}{
    \setboolean{normalnumbering}{false}
    \setboolean{sectnumbering}{true}}
\DeclareOption{noheader}{\setboolean{headers}{false}}
\DeclareOption{nocredit}{\setboolean{credit}{false}}
\DeclareOption{default}{\setboolean{fancyfont}{false}}
%------------------
\DeclareOption*{\ClassWarning{seto}{Unknown ‘\CurrentOption’}}
\ProcessOptions\relax
%---------------------------------------------
%colors, for mdframed
% Theme: big sur
% let's use mdframed already, its capabilities are a superset of those of tcolorbox
\definecolor{white}{HTML}{FFFFFF}% do not remove
\definecolor{black}{HTML}{000000}% do not remove
\definecolor{nyandarkgrn}{HTML}{0D534F}
\definecolor{probcolor}{HTML}{E5A29B} % brownish red?
\definecolor{theocolor}{HTML}{9BCAE8} % light blu
\definecolor{claimcolor}{HTML}{E6CE88} % yellow
\definecolor{remarkcolor}{HTML}{0D534F} % nyandarkgrn atm
\definecolor{quotecolor}{HTML}{91B9DC} % slightly darker blue
% (change hex codes to your taste to change box colors, its done rather automatically)
%-------------------------------------------------
%fonts, text formatting
\ifthenelse{\boolean{fancyfont}}{\RequirePackage[cmintegrals,cmbraces]{newtxmath}
	\RequirePackage{ebgaramond-maths}}{}
\DeclareRobustCommand{\sffamily}{\ifthenelse{\boolean{fancyfonts}}{\lato}{\cmss}}
\renewcommand{\emph}[1]{\color{nyandarkgrn}\sffamily\bfseries}
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
\newcommand{\exmp}[2][]{
    {\sffamily\bfseries Example}{\sffamily\optlabel{#1}}{\sffamily\bfseries.}\enskip#2\end{mdframed}}
\newcommand{\claim}[1]{
    \begin{mdframed}[style=claimdefault]{\sffamily\bfseries Claim ---\enskip}#1\end{mdframed}}
\newcommand{\lemma}[1]{
    \begin{mdframed}[style=claimdefault]{\sffamily\bfseries Lemma ---\enskip}#1\end{mdframed}}
\newcounter{theonum}
\newcommand{\thmlabel}{\ifthenelse{\boolean{normalnumbering}}{\thetheonum}{\thesection.\thetheonum}}
\newcommand{\theorem}[1]{\stepcounter{theonum}\begin{mdframed}[style=theodefault]
    {\sffamily\bfseries Theorem \thmlabel.\enskip}#1\end{mdframed}}
\newcommand{\theo}[1]{\begin{mdframed}[style=theodefault]
    {\sffamily\bfseries Theorem.\enskip}#1\end{mdframed}}
\newcounter{defnum}
\newcommand{\deflabel}{\ifthenelse{\boolean{normalnumbering}}{\thedefnum}{\thesection.\thedefnum}}
\newcommand{\definition}[1]{\stepcounter{defnum}\begin{mdframed}[style=theodefault]
    {\sffamily\bfseries Definition \deflabel.\enskip}#1\end{mdframed}}
\newcommand{\defn}[1]{\begin{mdframed}[style=theodefault]
    {\sffamily\bfseries Definition.\enskip}#1\end{mdframed}}
\newcommand{\remark}[1]{\begin{mdframed}[style=remarkdefault]
    {\color{nyandarkgrn}\sffamily\itshape Remark.\enskip}#1
    \end{mdframed}}
\newcommand{\note}[1]{\begin{mdframed}[style=remarkdefault]
    {\color{nyandarkgrn}\sffamily\itshape Note.\enskip}#1
    \end{mdframed}}
\newcommand{\pro}[1]{\begin{proof}#1\end{proof}}
%---------------------------------------------
%math operators/commands, from nyanya.sty, slight additions made
%delimiters
\newcommand{\parenth}[1]{\left( #1\right)}
\newcommand{\sqbrack}[1]{\left[ #1\right]}
\newcommand{\braces}[1]{\left\{ #1\right\}}
\newcommand{\abs}[1]{\left\lvert #1\right\rvert}
\newcommand{\norm}[1]{\left\lVert #1\right\rVert}
%floor/ceil; im so smart idk what \left and \right are, smh
\newcommand{\floor}[1]{\left\lfloor #1\right\rfloor}
\newcommand{\ceil}[1]{\left\lceil #1\right\rceil}
% subsitute for {a \brack b}
\newcommand{\sqbinom}[2]{\genfrac[]{0pt}{0}{#1}{#2}}
% operators, first two from lucky.sty, misc otherwise
\DeclareMathOperator{\cycsum}{\sum_{\text{cyc}}}
\DeclareMathOperator{\symsum}{\sum_{\text{sym}}}
\DeclareMathOperator{\cycprod}{\prod{\text{cyc}}}
\DeclareMathOperator{\symprod}{\prod{\text{sym}}}
\DeclareMathOperator{\dir}{\text{dir}}
\renewcommand{\deg}{^\circ}
\newcommand{\pow}[2]{\text{Pow}_{#1}(#2)}
\newcommand{\evalat}[2][]{\biggr\rvert_{#2}^{#1}}
\newcommand{\solqed}{\noindent\hfill$\blacksquare$}
\newcommand{\hbeq}[1]{\overset{#1}{=}}
\newcommand{\dangle}{\measuredangle}
%sets
\newcommand{\N}{\mathbb{N}}
\newcommand{\Z}{\mathbb{Z}}
\newcommand{\Q}{\mathbb{Q}}
\newcommand{\R}{\mathbb{R}}
\newcommand{\C}{\mathbb{C}}
\renewcommand{\P}{\mathbb{P}}% for example, PR^2
\newcommand{\F}{\mathbb{F}}
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
\setcounter{tocdepth}{3}
\newcommand{\toc}{
{\hypersetup{
	citecolor=black,
	filecolor=nyandarkgrn,
	linkcolor=black,
	urlcolor=black}
\emph{\tableofcontents}\newpage}}
%---------------------------------------------
% title customization
\RequirePackage{titling,titlesec,setspace}
\ifthenelse{\boolean{credit}}{\apptocmd{\maketitle}{\begin{center}
    {\color{nyandarkgrn}\Alegreya\itshape Produced using seto.cls}
    \end{center}}{}{}}{}
\patchcmd{\@maketitle}{\@title}{{\AlegreyaMedium\color{nyandarkgrn}\@title}}{}{}
\patchcmd{\@maketitle}{\@author}{{\Alegreya\color{nyandarkgrn}\@author}}{}{}
\patchcmd{\@maketitle}{\@date}{{\Alegreya\@date}}{}{}
%---------------------------------------------
% regular (sub)(sub)sections
\titleformat{\section}{\AlegreyaMedium\Large\color{nyandarkgrn}}
	{{\S\thesection}}{1em}{}
\apptocmd{\@section}{\ifthenelse{\boolean{sectnumbering}}{
        \setcounter{probnum}{0}\setcounter{exmpnum}{0}
        \setcounter{theonum}{0}\setcounter{defnum}{0}
        \setcounter{equation}{0}}{}}{}{}
\titleformat{\subsection}{\sffamily\bfseries\large\color{nyandarkgrn}}{
    {\thesubsection}}{1em}{}
\titleformat{\subsubsection}{\sffamily\bfseries\color{nyandarkgrn}}{
    {\thesubsubsection}}{1em}{}
%---------------------------------------------
% misc
%-----------
% exercises
\newcounter{exenum}
\newcommand{\exelabel}{\ifthenelse{\boolean{sectnumbering}}{\thesection.\theexenum}{\theexenum}}
\newcommand{\exercise}[1][]{\stepcounter{exenum}
	\par\noindent\emph{Problem \exelabel}{\color{nyandarkgrn}\sffamily\optlabel{#1}}\emph{.}\enskip}
%-----------
% sol
\newcommand{\solution}[1][]{\par\noindent {\sffamily\itshape Solution\optlabel{#1}.\enskip}}
%----------
% hpc problem numbering on school hw, because im lazy af
\newcounter{hpcExe}
\newcommand{\hpcExe}[1]{\stepcounter{hpcExe}\par\emph{\thehpcExe/#1; }}
\newcounter{hpcTopic}
\newcommand{\hpcTopic}[1]{\stepcounter{hpcTopic}\textit{Topic \thehpcTopic, #1.}}
\newcommand{\hpcContinueTopic}{\textit{Topic \thehpcTopic, continued.}}
\newcounter{revSectnum}
\newcommand{\hpcRevsect}[1][]{\stepcounter{revSectnum}\textit{Section \#\therevSectnum\ifthenelse{\isempty{#1}}{}{ (#1)}.}}
%----------
% 'fakecmd'
\newcommand{\fakecmd}[1]{\texttt{\textbackslash #1}}
%---------------------------------------------
%quotes
\newenvironment{quoteBox}{\begin{mdframed}[
	innerleftmargin=0.4cm,innertopmargin=0.25cm,innerbottommargin=0.25cm,innerrightmargin=0.4cm,
	leftmargin=3.25in, skipabove=0cm,skipbelow=0cm,
    linewidth=2.5pt,linecolor=quotecolor,backgroundcolor=quotecolor!25!white,
    topline=false,bottomline=false,rightline=false,leftline=true}]}
	{\end{mdframed}}
\renewcommand{\quote}[2]{\begin{quoteBox} {\Alegreya #1}
    {\begin{flushright}{\color{nyandarkgrn}\Alegreya #2}\end{flushright}}
    \end{quoteBox}}
%---------------------------------------------
