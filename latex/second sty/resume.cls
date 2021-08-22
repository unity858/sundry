%---------------------------------------------
%info
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{resume}[2021/08/22]
%% package name: resume.cls
%% Copyright 2021 by Neal Yan
%
% "let's modify seto.cls and actually create a cv already"
%
% This program can redistributed and/or modified under the terms 
% of the LaTeX Project Public License Distributed from CTAN 
% archives in directory macros/latex/base/lppl.txt; either 
% version 1 of the License, or (at your option) any later version.
%
%---------------------------------------------
\LoadClass[noheader,nogeo]{seto}
\RequirePackage[margin=2cm]{geometry}
\RequirePackage{xstring,xcolor}
\definecolor{maincolor}{HTML}{000000}
\definecolor{cvdarkgrn}{HTML}{0B524B}
\definecolor{misc}{HTML}{0B524B}
%---------------------------------------------
% robust stuff
\DeclareRobustCommand{\cvdarkgrn}{\color{cvdarkgrn}}
\DeclareRobustCommand{\misc}{\color{misc}}
\DeclareRobustCommand{\maincol}{\color{maincolor}}
%--------------------------------------------
% sectioning
\RequirePackage{titling,titlesec,setspace}
\newcommand{\sctn}[1]{{\par\noindent\Large\cvdarkgrn\AlegreyaMedium #1
    \vspace{0.1cm}\hrule\vspace{0.2cm}}}
\newcommand{\subsctn}[1]{{\par\noindent\AlegreyaMedium\large\cvdarkgrn #1\par\noindent}}
\newcommand{\subsubsctn}[1]{{\par\noindent\AlegreyaMedium\cvdarkgrn #1\par\noindent}}
%---------------------------------------------
\RequirePackage{mdframed}
\newcommand{\name}[1]{\def\realname{#1}}
\newcommand{\hugeName}{\begin{center}
    {\Huge\cvdarkgrn\lato\bfseries\realname}
    \end{center}}
\global\mdfdefinestyle{quote}{
    leftmargin=5.5cm,rightmargin=5.5cm,linecolor=white,backgroundcolor=white}
\newcommand{\cvQuote}[2]{\begin{mdframed}[style=quote]{\misc\Alegreya\itshape #1}
    \begin{flushright}{\misc\Alegreya #2}\end{flushright}
    \end{mdframed}}
%---------------------------------------------
% personal info links
\RequirePackage{fontawesome}
\newcommand{\links}[4]{\begin{center}\faLink\enskip #1\quad
    \faEnvelope\enskip #2\quad
    \faLinkedin\enskip #3\quad
    \faGithub\enskip #4\end{center}
    } % in order, website, mail, linkedin, github
%---------------------------------------------
% cv item, description, reference
\newcommand{\cvItem}[2]{\vspace{0.05cm}\par\noindent
    {\lato\bfseries\cvdarkgrn #1}\hfill {\lato #2}}
\newcommand{\desc}[1]{\par\noindent {\maincolor #1}}
\newcommand{\cvRef}[1]{\par\noindent(Reference: #1)}
%---------------------------------------------
