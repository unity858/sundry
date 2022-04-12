%---------------------------------------------
%info
% Beamer class by me
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{xuan}[beamer]
\LoadClass[9pt]{beamer}
\usetheme{Warsaw}
%% Last updated 22 August 2021
%% package name: xuan.cls
%% Copyright 2021 by Neal Yan
%
% This program can redistributed and/or modified under the terms 
% of the LaTeX Project Public License Distributed from CTAN 
% archives in directory macros/latex/base/lppl.txt; either 
% version 1 of the License, or (at your option) any later version.
%-----------------------------
% setup
\RequirePackage[nogeo,nohyperref,noimakeidx]{universal}
\setbeamersize{text margin left=0.5in,text margin right=0.5in}
%-----------------------------
% colors
\RequirePackage{xstring,xcolor}
\definecolor{white}{HTML}{FFFFFF}% do not remove
\definecolor{black}{HTML}{000000}% do not remove
\definecolor{tropical}{HTML}{64DFD7}
\definecolor{theocolor}{HTML}{CDFEFE} % light blu
\definecolor{claimcolor}{HTML}{F0F09C} % yellowish
%-----------------------------
% colors
\setbeamercolor*{structure}{bg=tropical!20!white,fg=tropical}
\setbeamercolor*{palette primary}{use=structure,fg=tropical,bg=structure.fg}
\setbeamercolor*{palette secondary}{use=structure,fg=white,bg=yellow}
\setbeamercolor*{palette tertiary}{use=structure,fg=white,bg=structure.fg!50!black}
\setbeamercolor*{palette quaternary}{fg=white,bg=tropical!60!white}
% title box & text
\setbeamercolor*{titlelike}{parent=tropical!80!white,bg=tropical,fg=white}
\setbeamercolor{frametitle}{bg=gray!10!white,fg=tropical}
\setbeamercolor{frametitle right}{bg=gray!60!white}
%
\setbeamercolor{section in toc}{fg=black,bg=white}
% block color config
\setbeamercolor{block title}{bg=tropical,fg=white}
\setbeamercolor{block body}{bg=block title.bg!15!white}
\setbeamercolor{block title alerted}{bg=claimcolor,fg=white}
\setbeamercolor{block body alerted}{bg=block title alerted.bg!15!white}
\setbeamercolor{block title example}{bg=theocolor,fg=white}
\setbeamercolor{block body example}{bg=block title example.bg!15!white,fg=black}

%-----------------------------
% boxes
\RequirePackage[most]{tcolorbox}
\RequirePackage[framemethod=TikZ]{mdframed}
% \newenvironment{prob}{\begin{tcolorbox}[title=Problem,colframe=probcolor,arc=1mm,left=5pt,right=5pt,top=4pt,bottom=4pt,colback=probcolor!25!white]}{\end{tcolorbox}}
\newenvironment{pt}[1][]{\begin{tcolorbox}[title=#1,colframe=tropical,arc=1mm,left=5pt,right=5pt,top=3pt,bottom=3pt,colback=tropical!20!white]}{\end{tcolorbox}}
%-----------------------------
