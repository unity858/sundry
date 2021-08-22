%---------------------------------------------
%info
% Beamer class by me
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{xuan}[beamer]
\LoadClass{beamer}
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
\definecolor{setogrn}{HTML}{0D534F}
\definecolor{probcolor}{HTML}{E5A29B} % brownish red?
\definecolor{theocolor}{HTML}{9BCAE8} % light blu
\definecolor{claimcolor}{HTML}{E6CE88} % yellow
\definecolor{remarkcolor}{HTML}{0D534F} % setogrn atm
%-----------------------------
% boxes
\RequirePackage[most]{tcolorbox}
\RequirePackage[framemethod=TikZ]{mdframed}
\newenvironment{prob}{\begin{tcolorbox}[title=Problem,colframe=probcolor,arc=1mm,left=5pt,right=5pt,top=4pt,bottom=4pt,colback=probcolor!25!white]}{\end{tcolorbox}}
\newenvironment{point}[1][]{\begin{tcolorbox}[title=#1,colframe=setogrn,arc=1mm,left=5pt,right=5pt,top=3pt,bottom=3pt,colback=setogrn!20!white]}{\end{tcolorbox}}
%-----------------------------
