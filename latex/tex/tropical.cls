%---------------------------------------------
%info
% Beamer class by me
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{tropical}[beamer]
\LoadClass[9pt,compress]{beamer}
\usetheme{Warsaw}
%% Last updated 12 April 2022
%% package name: tropical.cls
%% Copyright 2022 by Neal Yan
%
% This program can redistributed and/or modified under the terms 
% of the LaTeX Project Public License Distributed from CTAN 
% archives in directory macros/latex/base/lppl.txt; either 
% version 1 of the License, or (at your option) any later version.
%----------------------------------------
% setup
\RequirePackage[nogeo,nohyperref,noimakeidx]{universal}
\setbeamersize{text margin left=1.3cm,text margin right=1.3cm}
\addtobeamertemplate{block begin}
    {\setlength{\textwidth}{1.01\textwidth}}{}
\addtobeamertemplate{block alerted begin}{
  \setlength{\textwidth}{\textwidth}}{}
\addtobeamertemplate{block example begin}{
  \setlength{\textwidth}{\textwidth}}{}
%----------------------------------------
% colors
\RequirePackage{xstring,xcolor}
\definecolor{white}{HTML}{FFFFFF}% do not remove
\definecolor{black}{HTML}{000000}% do not remove
\definecolor{tropical}{HTML}{64DFD7}
\definecolor{exmpcolor}{HTML}{A2EBBE} % lig
\definecolor{claimcolor}{HTML}{F0F09C} %
\definecolor{skyblue}{HTML}{CDFEFE}
\definecolor{nyanpink}{HTML}{F3D9F1}
\definecolor{theocolor}{HTML}{B7E0FF}
%----------------------------------------
% colors
\setbeamercolor*{structure}{bg=tropical!80!white,fg=tropical}
% controls rfoot
\setbeamercolor*{palette primary}{use=structure,bg=tropical,fg=white}
% controls idk what
\setbeamercolor*{palette secondary}{use=structure,fg=white,bg=tropical}
% controls?
\setbeamercolor*{palette tertiary}{use=structure,fg=tropical,bg=skyblue}
% controls lfoot
\setbeamercolor*{palette quaternary}{use=structure,fg=black,bg=skyblue}
% title box & text
\setbeamercolor*{titlelike}{parent=skyblue,bg=tropical,fg=white}
% lhead does smth nontrivial
\setbeamercolor{frametitle}{bg=skyblue,fg=black}
% rhead fg does nothing
\setbeamercolor{frametitle right}{bg=tropical,fg=black}
% ??
\setbeamercolor{section in toc}{fg=black,bg=black}
% colors of blocks; bg=fill, bg=text
\setbeamercolor{block title}{bg=tropical,fg=white}
\setbeamercolor{block body}{bg=block title.bg!25!white}
\setbeamercolor{block title alerted}{bg=nyanpink!95!black,fg=white}
\setbeamercolor{block body alerted}{bg=block title alerted.bg!30!white,fg=black}
\setbeamercolor{block title example}{bg=exmpcolor,fg=white}
\setbeamercolor{block body example}{bg=block title example.bg!25!white,fg=black}
%----------------------------------------
% outline
\newcommand\outline[1][]{\begin{frame}{Outline}#1\tableofcontents\end{frame}}
% slide numbering
\newcommand*\oldmacro{}%
\let\oldmacro\insertshorttitle%
\renewcommand*\insertshorttitle{%
   \oldmacro\hfill%
   \insertframenumber\,/\,\inserttotalframenumber}
%----------------------------------------
% boxes
\RequirePackage[skins,most]{tcolorbox}
\RequirePackage[framemethod=TikZ]{mdframed}
% \newenvironment{prob}{\begin{tcolorbox}[title=Problem,colframe=probcolor,arc=1mm,left=5pt,right=5pt,top=4pt,bottom=4pt,colback=probcolor!25!white]}{\end{tcolorbox}}
\tcbuselibrary{breakable}
\tcbuselibrary{hooks}
\usetikzlibrary{shadings}
\newenvironment{claimBox}{\begin{mdframed}[
    innerleftmargin=0.3cm,innertopmargin=0.2cm,innerbottommargin=0.2cm,innerrightmargin=0.3cm,
    linecolor=claimcolor,linewidth=3pt,backgroundcolor=claimcolor!25!white,
    topline=false,bottomline=false,rightline=false,leftline=true]}{\end{mdframed}}
\newenvironment{thmBox}[1][]{\begin{tcolorbox}
    [enhanced,oversize,
    arc=1.5mm,left=1pt,right=1pt,top=1pt,bottom=1pt,toprule=0pt,
    toptitle=0pt,bottomtitle=0pt,titlerule=0pt,
    colframe=theocolor!95!black,colback=theocolor!30!white,
    bottomrule=0mm,leftrule=0mm,rightrule=0mm,
    drop fuzzy shadow=black,adjusted title=#1]}
    {\end{tcolorbox}}
\newcommand\claim[2][]
    {\begin{claimBox}{\bfseries Claim{\normalfont\optlabel{#1}}:}
    \enskip#2\end{claimBox}}
\newenvironment{thmblock}[1][Theorem]
    {\begin{thmBox}[#1]}{\end{thmBox}}
% env abbreviations
\newenvironment{thm}[1][]{\begin{thmblock}[#1]}{\end{thmblock}}
\newenvironment{warn}[1][]{\begin{alertblock}{#1}}{\end{alertblock}}
\newenvironment{exmp}[1][]{\begin{exampleblock}{#1}}{\end{exampleblock}}
%-----------------------------
