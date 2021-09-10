%---------------------------------------------
%info
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{amts}[2021/08/07]
\LoadClass[noheader,nogeo]{seto}
%% package name: usamts.cls 
%% Copyright 2021 by Neal Yan
%
% mock usamts cls, basically a modified version of seto.cls that can be used on usamts
%
% This program can redistributed and/or modified under the terms 
% of the LaTeX Project Public License Distributed from CTAN 
% archives in directory macros/latex/base/lppl.txt; either 
% version 1 of the License, or (at your option) any later version.
%
%---------------------------------------------
% amts info, geometry
\RequirePackage[left=1in,right=1in,top=1.2in,bottom=1in]{geometry}
\newcommand{\realname}[1]{\def\rname{#1}}
\newcommand{\username}[1]{\def\uname{#1}}
\newcommand{\amtsid}[1]{\def\id{#1}}
\newcommand{\amtsyr}[1]{\def\yr{#1}}
\newcommand{\amtsround}[1]{\def\round{#1}}
\newcounter{amtsprob}
%---------------------------------------------
% page setup, modified version of official;
% fonts/colors changed to the "seto" style;
\RequirePackage{fancyhdr}
\pagestyle{fancy}
\fancyhf{}
\setlength{\headheight}{0.75in}
\setlength{\oddsidemargin}{0in}
\setlength{\evensidemargin}{0in}
\setlength{\voffset}{-.5in}
\setlength{\headsep}{10pt}
\setlength{\textwidth}{6.5in}
\setlength{\headwidth}{6.5in}
\setlength{\textheight}{9in}
\fancyhead[L]{
    {\small\setogrn\Alegreya Student Name:\enskip\rname}\\
    {\small\setogrn\Alegreya Username:\enskip\uname}\\
    {\small\setogrn\Alegreya Student ID:\enskip\id}}
\fancyhead[C]{\Large\setogrn\Alegreya USAMTS Solutions}
\fancyhead[R]{\begin{tabular}{c|c|c}
     {\small\emph{Year}} & {\small\emph{Round}} & {\small\emph{Problem}} \\\hline
     {\small\Alegreya\yr} & {\small\Alegreya\round} & {\small\Alegreya\theamtsprob}
\end{tabular}}
\fancyfoot[L]{}
\fancyfoot[C]{}
\fancyfoot[R]{{\setogrn\Alegreya Page \thepage\ of Problem \theamtsprob}}
\renewcommand{\headrulewidth}{0.5pt}
\renewcommand{\footrulewidth}{0.5pt}
\setlength{\textwidth}{6.5in}
%---------------------------------------------
% amtssol env
\newenvironment{amtssol}{
    \stepcounter{amtsprob}
    \setcounter{page}{1}
    \setcounter{equation}{0}
    \setcounter{probnum}{0}
    \setcounter{exmpnum}{0}
    \setcounter{theonum}{0}
    \setcounter{defnum}{0}
    \section{Problem \theamtsprob}}
    {\eject}
%---------------------------------------------
