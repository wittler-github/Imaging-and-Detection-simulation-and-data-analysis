\ProvidesClass{my_report}
\LoadClass[11pt, titlepage]{article}
\usepackage[paperwidth=21cm,paperheight=29.7cm,includehead,headheight=2.0cm,pdftex,hmargin={3cm,2.5cm},vmargin={0cm,2cm},]{geometry} 
\setlength{\parindent}{0cm} 
\renewcommand{\baselinestretch}{1.2}
\parskip=1em

\renewcommand{\contentsname}{Table of Contents}

\renewcommand{\normalsize}{\fontsize{11pt}{11pt}\selectfont}
\renewcommand{\title}[1]{\textsf{\Huge{#1}}}
\renewcommand{\author}[1]{\textsf{\Huge{#1}}}
\newcommand{\subtitle}[1]{\textsf{\textbf{\Large{#1}}}}
\newcommand{\abstractheading}[1]{\textsf{\textbf{\LARGE{#1}}}}
\newcommand{\code}[1]{\texttt{\footnotesize{#1}}}

\newcommand{\codename}{NADIA.0.5.tar}
\newcommand{\name}{{\tt DIXE~}}
\newcommand{\ver}{0.4}
\newcommand{\website}{http://coecxs.ph.unimelb.edu.au/}

%\makeatletter
%\renewcommand*\l@subsection{\@dottedtocline{2}{1.5em}{2.0em}}
%\renewcommand*\l@subsubsection{\@dottedtocline{3}{3.5em}{3.0em}}
%\makeatother

\RequirePackage[dvips]{color}
\definecolor{mycolour}{rgb}{0.05,0.2,0.5}
\definecolor{gray}{rgb}{0.5,0.5,0.5}

\renewcommand{\section}{\@startsection
{section}%                    % the name
{1}%                          % the level
{0mm}%                        % the indent
{10.6mm}%             % the beforeskip
{4.2mm}%           % the afterskip
{\LARGE\bfseries\sffamily\color{mycolour}}}  % the style

\renewcommand{\subsection}{\@startsection
{subsection}%                    % the name
{2}%                          % the level
{0mm}%                        % the indent
{6.4mm}%             % the beforeskip
{1.1mm}%           % the afterskip
{\Large\bfseries\sffamily\color{mycolour}}}  % the style

\renewcommand{\subsubsection}{\@startsection
{subsubsection}%                    % the name
{3}%                          % the level
{0mm}%                        % the indent
{4.2mm}%             % the beforeskip
{1.1mm}%           % the afterskip
{\normalsize\bfseries\sffamily\color{mycolour}}}  % the style

%% Chapter headings should be centered, uppercase, and at the top of the page.
\def\@makechapterhead#1{%
  { \parindent \z@ \raggedright \normalfont
    %\centering
    \ifnum \c@secnumdepth >\m@ne
        \huge\textsf{\@chapapp\space \thechapter:}
        % \par\nobreak
        %\vskip 20\p@
    \fi
    \interlinepenalty\@M
    \huge \bfseries \textsf{#1}\par\nobreak
    \rule{5cm}{0.5pt}
    \vskip 20\p@
  } }
  
\def\@makeschapterhead#1{%
  %\vspace*{50\p@}%
  { \parindent \z@ \raggedright
    %\centering
    \normalfont
    \interlinepenalty\@M
    \huge \bfseries  \textsf{#1}\par\nobreak
    \rule{5cm}{0.5pt}
    \vskip 20\p@

  }}
  
 \renewenvironment{abstract}{%
      \chapter*{\abstractname}%
      \addcontentsline{toc}{chapter}{\abstractname}
 }
     
     
\makeatletter
\renewcommand{\l@section}{\bfseries\@dottedtocline{1}{0em}{2.3em}}
\renewcommand{\l@subsection}{\normalfont\@dottedtocline{2}{2em}{2.3em}}
\renewcommand{\l@subsubsection}{\normalfont\@dottedtocline{3}{3.5em}{2.3em}}
%\renewcommand{\l@subsubsection}{\normalfont\@dottedtocline{4}{2em}{2.3em}}
\makeatother

\def\maketitle{\begin{titlepage}
\thispagestyle{empty}
\let\footnotesize\small \let\footnoterule\relax \setcounter{page}{0}


\begin{center}
\vskip 90 mm

\includegraphics[height=3cm]{cxs_logo.eps}\\[1cm]

%\rule{4cm}{1pt}

\title{\color{mycolour}{\textbf{DIXE: Diffractive Imaging for X-rays and Electrons}\\[0.5cm] Software Manual\\[0.5cm] Version \ver \\[0.5cm] \rule{4cm}{1pt}\\[0.7cm]}
}

\author{\Large N. Davidson$^1$, B. Abbey$^1$, H. Quiney$^1$, M. Jones$^2$, \\[0.5cm] 
C. Putkunz$^1$, A. Torrance$^1$, H. Wittler$^2$, \\[0.5cm] 
A. Carroll$^2$, M. Luu$^2$  and K. Nugent$^1$}\\[1.0cm]

{\em $^1$ University of Melbourne, Victoria, Australia.  }\\[0.5cm]
{\em $^2$ La Trobe University, Victoria, Australia.  }\\


\vskip 15 mm

\textbf\sf{cxs-software@physics.unimelb.edu.au} 

\end{center}
\vskip 10 mm

\centerline{\normalsize \today}


\end{titlepage}}

\pagestyle{plain}

