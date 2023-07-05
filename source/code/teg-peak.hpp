#ifndef TEGPEAK_HPP
#define TEGPEAK_HPP

#include <stdlib.h>     // For use of exit command
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstring>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/encoding.h>
#include <libxml/xmlwriter.h>
#include <libxml/xpathInternals.h>
#include "random.hpp"
using namespace std;

#define INFIN 2147483647

/* Code for the TEGP generator. All documentation are written using LaTeX in the
   sources and can be generated using doxygen.

   Author: Lars Relund Nielsen - lars@relund.dk
   Version 1.66
*/



/** \mainpage Description

\latexonly

\title{The TEGP generator}

\author{%
Lars Relund Nielsen\thanks{Corresponding author (\url{lars@relund.dk}).}\\
{\small CORAL - Cluster for OR Applications in Logistics, Department of Economics and Business,
School of Business and Social Sciences, Aarhus University, Fuglesangs Allé 4,
DK-8210 Aarhus V, Denmark.}
}%

\date{Version 1.66 - July 2013}

\maketitle

%% --------------------------------------------------------------------- %%
%% Abstract and Keywords
\begin{abstract}
This manual provides documentation of the \emph{Time-Expanded Generator with
Peaks} (\emph{TEGP}) which generates instances of stochastic time-dependent
networks. The program includes several features inspired by typical aspects of
road networks (congestion effects, waiting, random perturbations etc.).

\vskip 0.3cm

{\em Keywords:} stochastic time-dependent networks, testing, generator.
\end{abstract}
%% --------------------------------------------------------------------- %%

\pagenumbering{roman}

\tableofcontents

\pagenumbering{arabic}

\section{Introduction}

The \emph{Time-Expanded Generator with Peaks} (\emph{TEGP}) generates
instances of \emph{stochastic time-dependent networks} (\emph{STD networks}).
The program includes several features inspired by typical aspects of road
networks (congestion effects, waiting, random perturbations etc.). Note that
the TEGP generator like other generators only models a fraction of a real
network. However, it provides alternative choices that may affect the behavior
of the algorithms. For an short introduction to STD networks and the notation
used see \Fref{apdx:prelim}.

\begin{figure}[t]
   \begin{center}
   \includegraphics{graphic/grid_3x3}
   \end{center}
   \caption{A topological grid network ($3\times3$ grid).}
   \label{fig:grid}
\end{figure}

A topological grid network $G$ of \index{base}\emph{base }$b$ and
\index{height}\emph{height }$h$ is assumed, and we search for optimal
strategies from the bottom-right corner node (origin\emph{ o}) to the upper
left corner node (destination $d)$. This choice is motivated by the fact that
each origin-destination path has at least $b+h-2$ arcs, and the number of such
paths grows exponential with the size of $G$. A topological grid network $G$
of base $3$ and height $3$ is shown in \Fref[plain]{fig:grid}. Note that no
arcs enter the origin and no arcs exit the destination.

The generator considers \emph{cyclic} time periods (e.g.\ a day) and the time
horizon $H$ is the (finite) number of time instances in a cycle multiplied by
the number of cycles. In each cyclic period there are some
\index{peak}\emph{peak periods} (e.g.\ rush hours). Each peak consists of
three parts; a \emph{transient} part  where the traffic increases, a
\emph{pure peak} part  where the traffic stays the same and a transient part
where the traffic decreases again. Peaks are placed at the same time in each
cycle. If no peaks are wanted then the length of the peak may be set to zero.

%% ----------------------------------------------------------------------------
\section{Generating travel time distributions}
\label{sec:travel_time}

The travel time distributions are found by first generating the
\emph{off-peak} mean travel times $\mu(u,v)\in\{lb_{T},...,ub_{T}\}$ for all
arcs $(u,v)$.

Let $\mu_{uv}(t)$ denote the mean travel time for each possible leaving time
$t\in L(u,v)$. The mean $\mu_{uv}\left(  t\right)  $ follows a pattern like
the line in \Fref{fig:tegp_mean} where the mean for each leaving time is shown
with a circle. In off-peaks $\mu_{uv}\left( t\right) =\mu\left( u,v\right)$,
i.e. the constant mean travel time in off-peaks. At the beginning of a peak,
$\mu_{uv}\left( t\right) $ increases from $\mu\left( u,v\right)  $ to
$\mu\left(  u,v\right) \left( 1+\psi\right) $, where $\psi$ denotes the
\index{peak increase parameter}\emph{peak increase parameter}, then stays the
same during the pure peak period, and then decreases to $\mu\left(  u,v\right)
$ again. Which arc the peak effect is applied to is controlled by the flag $flag_{T}$
which may take the following values
\begin{enumerate}
  \item[0 -] The peak effect is applied to all arcs.
  \item[1 -] The peak effect is applied to north and west arcs only. The other arcs (south and east)
    use $\mu_{uv}\left( t\right) =\mu\left( u,v\right)$ for all leaving times.
  \item[2 -] The peak effect is applied to west and east arcs only. The other arcs
    use $\mu_{uv}\left( t\right) =\mu\left( u,v\right)$ for all leaving times.
\end{enumerate}
Option 1 and 2 have been provided to model roads which not are so ``Peak sensitive''.

\begin{figure}[t]
   \begin{center}
   \includegraphics[width=0.8\linewidth]{graphic/mean_travel_time}
   \caption{Mean travel time for an arc ($\psi=100\%$).}
   \end{center}
   \label{fig:tegp_mean}
\end{figure}

Let $X\left(  u,v,t\right)  $ denote the travel time to node $v$ when leaving
node $u$ at time $t$ along arc $\left(  u,v\right)$. Given the mean travel
time $\mu_{uv}(t)$ the travel time distribution is found as follows.

\begin{enumerate}
\item Find $\sigma_{uv}\left( t\right) =\rho\mu_{uv}\left(  t\right)  $ where $\rho$ is
the \index{standard deviation mean ratio}\emph{standard deviation mean ratio}.

\item The set of possible travel times is then $\{t_1,...,t_{\kappa(u,v,t)}\}=
\{\lfloor\mu_{uv}(t)-\sigma_{uv}(t)\rfloor,...,\mu_{uv}(t),...,\lceil\mu_{uv}(t)+\sigma_{uv}(t)\rceil\}$.
Moreover, only positive travel times in $\{t_1,...,t_{\kappa(u,v,t)}\}$ are
considered.

\item Note $X(u,v,t)=t_1+Y(u,v,t)$ where $Y$ is a discrete random variable
taking the values $\{0,...,\kappa(u,v,t)-1\}$. Given the positive travel times
$\{t_1,...,t_{\kappa(u,v,t)}\}$ we now set $Pr(X=t_i)=Pr(Y=i-1)$ where $Y\sim
bi(q-1,1/2)$. As a result $E(X)=\mu_{uv}(t)$.
\end{enumerate}

Note, using the setting above gives higher mean travel time and higher
standard deviation in peaks. Moreover, only the interval of possible off-peak
travel times $\left\{ lb_{T},...,ub_{T}\right\}  $ is given as an input
parameter to the TEGP generator. The actual interval of possible travel times
$I_{T}$ depends on the peak increase parameter $\psi$ and the standard
deviation mean ratio $\rho$, i.e.

\begin{equation} I_{T}=\left\{  \lfloor\left(  1-\rho\right)  lb_{T}\rfloor,\dots,\lceil\left(
1+\psi\right)  \left( 1+\rho\right)  ub_{T}\rceil\right\}  \label{eq:I_T}
\end{equation}

A flag $flag_{static}$ have been provided such that it is possible to use static travel time and
costs in off-peaks. If $flag_{static}=0$ travel times are generated as described above. If
$flag_{static}=1$ a deterministic travel time $\mu_{uv}\left( t\right) =\mu\left( u,v\right)$ is
used in off-peaks. This feature is provided if you only want to have a stochastic nature in peaks.
For costs the random pertubation is ignored in off-peaks if $flag_{static}=1$ (see \Fref{sec:travel_costs})

Finally, symmetric mean travel times can be used for the arcs in $G$. In this
case $\mu_{uv}\left( t\right)  =\mu_{vu}\left(  t\right)  $. Note symmetric mean travel times is not possible if $flag_T=1$

%% ----------------------------------------------------------------------------
\section{Generating costs}

Three different types of costs are considered, namely, travel costs, waiting
costs and penalty costs. Currently, negative costs are not accepted.

\subsection{Generating travel costs}
\label{sec:travel_costs}

Two costs $c_{i}\left(  u,v,t\right)$, $i=1,2$ for each arc $\left(
u,v\right)$ and leaving time $t\in H$ are generated, since we may need two
costs if bicriterion route choice is considered. The way the costs are
generated is controlled using two flags. The first flag, $flag_{cor}$, specify
the correlation between the two off-peak costs $c_{i}\left( u,v\right)$,
$i=1,2$. The second flag, $flag_C$, specify how the costs, given an arc, for
different leaving times depend on each other.

The following values of $flag_{cor}$ are possible:

\begin{enumerate}
  \item[0 -] both costs $c_{i}(u,v)$ are random in $\left\{
lb_{C},...,ub_{C}\right\} $.
  \item[1 -] $c_{2}(u,v) =  ub_C - (c_{1}(u,v)-lb_C)$
  \item[2 -] The costs are generated as follows:
  \begin{align*}
c_{1}\left(  u,v\right)   & <\frac{ub_{C}-lb_{C}}{2}\Rightarrow c_{2}\left(
u,v\right)  \in\left\{
ub_{C}-\left(  c_{1}\left(  u,v\right)  -lb_{C}\right) ,\dots,ub_{C}\right\} \\
c_{1}\left(  u,v\right)   &  \geq\frac{ub_{C}-lb_{C}}{2}\Rightarrow
c_{2}\left(  u,v\right)  \in\left\{  lb_{C},\dots,lb_{C}+\left(
ub_{C}-c_{1}\left( u,v\right)  \right)  \right\}
\end{align*}
  \item[3 -] NETMAKER costs see \citet{Skriver99b}.
          ``...if one cost is between 1 and 33,
           the other is between 67 and 100''
\end{enumerate}

Note for $flag_{cor}$ equal 1, 2 and 3 the two costs are negatively
correlated. This is a typical situation in hazardous material transportation,
where travel cost and risk/exposure are conflicting.

Off-peak costs generated in the interval $\left[  1,1000\right] \times\left[
1,1000\right] $ for the four correlation types are shown in \Fref
[plain]{fig:cor}.

\begin{figure}[tb]
   \centering
   \subfigure[$flag_{cor}=0$]{\includegraphics[width=0.45\linewidth]{graphic/off-peak_costs_cor0}}
   \subfigure[$flag_{cor}=1$]{\includegraphics[width=0.45\linewidth]{graphic/off-peak_costs_cor1}}\\
   \subfigure[$flag_{cor}=2$]{\includegraphics[width=0.45\linewidth]{graphic/off-peak_costs_cor2}}
   \subfigure[$flag_{cor}=3$]{\includegraphics[width=0.45\linewidth]{graphic/off-peak_costs_cor3}}
   \caption{Off-peak costs for the different correlation types.}
   \label{fig:cor}
\end{figure}

Given costs $c_{i}(u,v)$ the generation of costs $c_{i}\left(  u,v,t\right)$,
$i=1,2$ may take three components into account: the off-peak cost, the peak
effect, and a random perturbation. The \index{random perturbation}\emph{random
perturbation} introduces small variations in the cost, due to other factors
not intercepted by the peak, e.g.\ special information about the cost at
exactly that leaving time. Given a cost $c$, we generate a perturbation
$\xi\in\left[  -r_{\xi},r_{\xi}\right] $, where \index{random
perturbation!range of}range \emph{r}$_{\xi}$ is a small percentage. Then, the
cost after applying the random perturbation becomes $c(1+\xi)$. The following
values of $flag_{C}$ are possible:

\begin{enumerate}
\item[0 -] Costs are generated independently for each leaving time and
correlated as specified by $flag_{cor}$. That is, we have costs in the
interval $\{lb_{C},...,ub_{C}\}$. The random perturbation is not applied. A
plot of the costs for a specific arc is shown in \Fref{fig:costs_f0_cor2}. In
the left figure costs $c_{i}\left(  u,v,t\right)$, $i=1,2$ are displayed for
each leaving time (1. cost is the solid line, 2. cost is the dotted line) and
in the right figure the costs $(c_{1},c_2)$ are shown. Note that by using this
setting, e.g.\ the costs of leaving $u$ at time $t$ may be 10 and the costs of
leaving $u$ at time $t+1$ may be 500. In a road network this is probably not
realistic and peak dependent cost can be generated instead. However, random
costs can be used to see if e.g. a priori algorithms are robust, see
\citet{Relund04b}.

\item[1 -] Consider node $u$ having a west, north, east and south arc.\footnote{Some arcs are not considered in output if the arc do not exist in the underlying grid.}
First the off-peak costs $c_{i}(u,v) \in\{lb_{C},...,ub_{C}\}$, $i=1,2$, are
generated for the west and north arcs with a correlation as specified by
$flag_{cor}$. Let $c^{\min}_i$ denote the minimum cost generated for criterion
$i$ of the north and west arcs. For the east and south arc off-peak costs
$c_{i}(u,v) \in\{lb_{C},...,c_{i}^{\min}\}$, $i=1,2$, are now generated. In
this way, east and south arcs are expected to have smaller costs than the
corresponding west or north arcs. Finally, the costs $c_{i}\left(
u,v,t\right)$, $i=1,2$ are found by applying the random perturbation. Note can
only be used with non symmetric arcs. Costs of a west arc $(5,2)$ and an east
arc $(5,8)$ (of the grid in \Fref{fig:grid}) are given in
\Fref{fig:costs_f1_cor2}. Note the only difference in e.g. $c_1(5,2,t)$,
$t=1,...,H$ is due to the random perturbation (equal 10\%).

\item[2 - ] First the off-peak costs $c_{i}(u,v) \in\{lb_{C},...,ub_{C}\}$,
$i=1,2$, are generated with a correlation as specified by $flag_{cor}$. Next,
the costs $c_{i}\left(  u,v,t\right)$, $i=1,2$ are found by applying the
random perturbation. A plot of the costs for a specific arc are given in
\Fref{fig:costs_f2_cor2}.

\item[3 -] Generation of \emph{peak dependent costs}:
First, the off-peaks cost $c_{i}\left(  u,v\right) \in\left\{
lb_{C},...,ub_{C}\right\} $, $i=1,2$, are generated with a correlation as
specified by $flag_{cor}$. After generating the off-peak costs the peak effect
is taken into account. For each arc $\left( u,v\right)  $, the costs, if
leaving node $u$ at an off-peak time, are $\hat{c}_{i}\left(  u,v,t\right)
=c_{i}\left( u,v\right) ,$ $i=1,2$. At the beginning of a peak, the costs
$\hat{c}_{i}\left( u,v,t\right) $ increase from $c_{i}\left(  u,v\right)  $ to
$c_{i}\left( u,v\right)  \left( 1+\psi\right)  $, then stays the same during
the pure peak period, and then decrease to $c_{i}\left(  u,v\right)$ again.
Finally, the costs $c_{i}\left(  u,v,t\right)$, $i=1,2$ are found by applying
the random perturbation to $\hat{c}_{i}\left(  u,v,t\right)$. A plot of the
costs for a specific arc are given in \Fref{fig:costs_f3_cor2}.

\item[4 -] Peak dependent costs: The first cost act like for $flag_{C}=3$.
But the second cost decrease in peaks instead of increasing. A plot of the
costs for a specific arc is given in \Fref{fig:costs_f4_cor2}.

\item[5 -] Peak dependent costs: The off-peak costs $c_{i}\left(
u,v\right)$ are generated as under $flag_{C}=3$. If $c_{1}\left(
u,v\right)<c_{2}\left(  u,v\right)$ then a peak cost $c_{1}^p(u,v)$ is
generated randomly in $\{c_{1}\left( u,v\right),\ldots,c_{2}\left(
u,v\right)\}$ and a peak cost $c_{2}^p(u,v,)$ is generated randomly in
$\{c_{2}\left( u,v\right),\ldots,ub_C\}$. As a result $c_{1}\left(
u,v\right)\leq c_{1}^p(u,v,)\leq c_{2}\left( u,v\right)\leq c_{2}^p(u,v)$.
Similar is the case $c_{1}\left( u,v\right)\geq c_{2}\left( u,v\right)$. Next,
the costs $c_{i}\left( u,v,t\right)$, $i=1,2$ are found by applying the random
perturbation to $c_{i}(u,v)$ in off-peaks and $c_{i}^p(u,v)$ in peaks. A plot
of the costs for a specific arc is given in \Fref{fig:costs_f5_cor2}.

\item[6 -] Peak dependent costs as if $flag_{C}=3$ for horizontal arcs
and costs are generated as for $flag_{C}=2$ for vertical arcs. A plot of the
costs for a vertical and a horizontal arc  are given in
\Fref{fig:costs_f6_cor2}.
\end{enumerate}

\begin{figure}[ptb]
   \centering
   \includegraphics[width=\linewidth]{graphic/costs_flag0_cor2}
   \caption{Travel costs for $flag_C=0$ and $flag_{cor}=2$.}
   \label{fig:costs_f0_cor2}
\end{figure}

\begin{figure}[ptb]
   \centering
   \includegraphics[width=\linewidth]{graphic/costs_flag1_cor2}
   \caption{Travel costs for $flag_C=1$ and $flag_{cor}=2$.}
   \label{fig:costs_f1_cor2}
\end{figure}

\begin{figure}[ptb]
   \centering
   \includegraphics[width=\linewidth]{graphic/costs_flag2_cor2}
   \caption{Travel costs for $flag_C=2$ and $flag_{cor}=2$.}
   \label{fig:costs_f2_cor2}
\end{figure}

\begin{figure}[ptb]
   \centering
   \includegraphics[width=\linewidth]{graphic/costs_flag3_cor2}
   \caption{Travel costs for $flag_C=3$ and $flag_{cor}=2$.}
   \label{fig:costs_f3_cor2}
\end{figure}

\begin{figure}[ptb]
   \centering
   \includegraphics[width=\linewidth]{graphic/costs_flag4_cor2}
   \caption{Travel costs for $flag_C=4$ and $flag_{cor}=2$.}
   \label{fig:costs_f4_cor2}
\end{figure}

\begin{figure}[ptb]
   \centering
   \includegraphics[width=\linewidth]{graphic/costs_flag5_cor2}
   \caption{Travel costs for $flag_C=5$ and $flag_{cor}=2$.}
   \label{fig:costs_f5_cor2}
\end{figure}

\begin{figure}[ptb]
   \centering
   \includegraphics[width=\linewidth]{graphic/costs_flag6_cor2}
   \caption{Travel costs for $flag_C=6$ and $flag_{cor}=2$.}
   \label{fig:costs_f6_cor2}
\end{figure}

Symmetric travel costs may be used. However, note the costs are symmetric
before the random perturbation is applied, i.e.\ $\hat{c}_{i}\left(
u,v,t\right) =\hat{c}_{i}\left( v,u,t\right)  ,$ $i=1,2$ and then the random
perturbation is applied to obtain $c_{i}\left(  u,v,t\right)$. As a result we
do not have exact symmetry, costs may vary a bit. If you want exact symmetry
you may set the range $r_{\xi}$ of the random perturbation to zero. Symmetric
travel costs are controlled with the $flag_{sym}$ flag. Possible options are:

\begin{description}
  \item[0 -] Do not consider symmetry.
  \item[1 -] Consider symmetry.
\end{description}

\subsection{Generating waiting costs}
\label{sec:wait}

Waiting are considered if the input parameter $ub_{W}>0$. In this case waiting
costs $c_{i}\left( u,t,t+1\right) ,$ $i=1,2$ for each node $u$ in the grid
except for the origin $o$ and destination $d$ are generated.

How waiting costs are generated are specified by the $flag_W$ flag. We have 2
possible values:

\begin{description}
\item[0 - ] Node waiting cost $c_{i}\left( u\right)$ is generated in
$\{lb_{W},...,ub_{W}\}$ and correlated as specified by the $flag_{cor}$.
Waiting costs $c_{i}\left( u,t,t+1\right)$ are then generated by applying the
random perturbation. A plot of the waiting costs is given in
\Fref{fig:waiting_f0_cor2}.

\item[1 - ] Waiting costs $c_{i}\left( u,t,t+1\right)$ are generated independently
for each leaving time in $\{lb_{W},...,ub_{W}\}$ and correlated as specified
by the $flag_{cor}$. The random perturbation is not applied. A plot of the
waiting costs is given in \Fref{fig:waiting_f1_cor2}.
\end{description}

\begin{figure}[ptb]
   \centering
   \includegraphics[width=\linewidth]{graphic/waiting_flag0_cor2}
   \caption{Waiting costs for $flag_W=0$ and $flag_{cor}=2$.}
   \label{fig:waiting_f0_cor2}
\end{figure}

\begin{figure}[ptb]
   \centering
   \includegraphics[width=\linewidth]{graphic/waiting_flag1_cor2}
   \caption{Waiting costs for $flag_W=1$ and $flag_{cor}=2$.}
   \label{fig:waiting_f1_cor2}
\end{figure}

Note the peak effect is not used for waiting costs.


%% ----------------------------------------------------------------------------
\subsection{Generating penalty costs}
\label{sec:penalty}

Penalty costs for the destination node $d$ are generated if the input
parameter $ub_P>0$. Otherwise, the penalty costs are zero. How penalty costs
are generated is specified by the $flag_P$ flag. We have 4 possible values:

\begin{description}
\item[0 - ] Penalty costs are generated independently
for each arrival time in the interval $\{lb_{P},...,ub_{P}\}$ and correlated
as specified by the $flag_{cor}$. A plot of the penalty costs is given in
\Fref{fig:penalty_f0_cor2}.

\item[1 - ] Both costs penalizes arrivals in the middle of the time horizon
$H$. A plot of the penalty costs is given in \Fref{fig:penalty_f1} (first and
second cost equals the solid line).

\item[2 -] Both costs penalizes early/late arrivals.
A plot of the penalty costs is given in \Fref{fig:penalty_f2}.

\item[3 -] The second cost penalizes arrivals in the middle of the time horizon
the first cost penalizes early/late arrivals. A plot of the penalty costs is
given in \Fref{fig:penalty_f3}.
\end{description}

Note the random perturbation is NOT used when generating penalty costs and for
$flag_P>0$ the $flag_{cor}$ have no effect on the penalty costs.

\begin{figure}[ptb]
   \centering
   \includegraphics[width=\linewidth]{graphic/penalty_flag0_cor2}
   \caption{Penalty costs for $flag_P=0$ and $flag_{cor}=2$.}
   \label{fig:penalty_f0_cor2}
\end{figure}

\begin{figure}[ptb]
   \centering
   \includegraphics[width=\linewidth]{graphic/penalty_flag1}
   \caption{Penalty costs for $flag_P=1$.}
   \label{fig:penalty_f1}
\end{figure}

\begin{figure}[ptb]
   \centering
   \includegraphics[width=\linewidth]{graphic/penalty_flag2}
   \caption{Penalty costs for $flag_P=2$.}
   \label{fig:penalty_f2}
\end{figure}

\begin{figure}[ptb]
   \centering
   \includegraphics[width=\linewidth]{graphic/penalty_flag3}
   \caption{Penalty costs for $flag_P=3$.}
   \label{fig:penalty_f3}
\end{figure}




%% ----------------------------------------------------------------------------
\section{Calculating the time horizon}\label{sec:horizon}

The time horizon $H$ for the STD network may be given as an input parameter to
the TEGP generator (if has a positive value, see \Fref{sec:input}).

Note that the time horizon depends on the size of $G$ and the
possible travel times generated for the arcs in $G$. If the input parameter specifying $H$
is not positive then the time horizon is found using a
preprocessing step:

First, note that due to \eqref{eq:I_T} an upper bound on a possible travel
time for an arc is $ub=(1+\psi)(1+\rho)ub_{T}$. As result $H_{ub}=(b+h)ub$ is
an upper bound on the travel time of a path of length $b+h$. Upper bound $ub$
is in general not very tight since there may be large fluctuations in travel
time for different arcs and leaving time. Hence we use an estimate for the
average maximal travel time for a path. This is done by generating all travel
time distributions for each arc and leaving time for time horizon $H_{ub}$.
Now by storing the maximum possible travel time for each distribution we can
calculate the average maximum possible travel time $ub_{ave}$. We now set the
time horizon H to
\[
H=(b+h)ub_{ave}
\]
As a result ``roughly'' all $o--d$ paths of length $b+h$ can be traveled in
the time horizon. Moreover, note that often there exists many paths containing
arcs with maximum possible travel time below $ub_{ave}$. Therefore paths of
length greater than $b+h$ may often be possible to travel in the time horizon.

%% ----------------------------------------------------------------------------
\section{Running the program}

The program uses command line passing for catching a number of run directives.
The following flags/options can be used:

\begin{description}
  \item[-verbose] print out a lot of information to standard output (optional).
  \item[-out] file name of the output file follows (without extension).
  \item[-xml] Output file generated in xml format (see \Fref{sec:output}).
  \item[-f5] Output file generated in time-expanded hypergraph format (f5 format - see \Fref{sec:output}).
\end{description}

Input is read from standard input which can be piped.

Example: \texttt{tegp < input -out output -xml}


%% ----------------------------------------------------------------------------
\section{Input parameters}\label{sec:input}

Input Parameters must be integers and are read from standard input or piped
using a file with parameters separated by spaces as illustrated
below\footnote{The parameters may also be given on a single line.}.

\begin{tabularx}{\linewidth}{lX}
\hline
   $b$ $h$            & Base and height in grid. \\
   $H_{cycle}$        & Time instances in a cycle. \\
   $p$                & Number of peaks in a cycle. \\
   $t_{trans}$        & Time instances in transient peak. \\
   $t_{pure}$         & Time instances in pure peak. \\
   $t_p$              & First peak starting time. \\
   $H$                & Time horizon used (if positive). If negative the time horizon is calculated by the generator (see \Fref{sec:horizon}). \\
   $\psi$             & Mean travel time/cost increase in peaks (pct). \\
   $\rho$             & Variance/mean ratio (pct). \\
   $lb_{P}$ $ub_{P}$  & Min and max penalty cost. If $ub_P<0$ then no penalty costs are generated. \\
   $flag_{P}$         & Dependencies of the penalty costs (see \Fref{sec:penalty}).  \\
   $lb_{T}$ $ub_T$    & The mean traveltime interval. \\
   $flag_T$           & Peak effect on arcs fag (see \Fref{sec:travel_time}). \\
   $lb_{W}$ $ub_{W}$  & Min and max waiting cost. If $ub_W<0$ then no waiting arcs are generated. \\
   $flag_W$           & Waiting costs flag (see \Fref{sec:wait}). \\
   $lb_{C}$ $ub_{C}$  & Min and max travel cost.  \\
   $flag_{cost}$      & Dependencies of the travel costs (see \Fref{sec:travel_costs}). \\
   $flag_{sym}$       & Travel time symmetry flag (see \Fref{sec:travel_costs}).  \\
   $flag_{cor}$       & Correlation type between costs (see \Fref{sec:travel_costs}).\\
   $rand$             & Random perturbation (promille, see \Fref{sec:travel_costs}). \\
   $seed$             & Seed. \\
\hline
\end{tabularx}

Note that peaks are distributed evenly in the time interval $[t_p;H_{cycle}]$.

%% ----------------------------------------------------------------------------
\section{Output}\label{sec:output}

\begin{figure}[tbp]
\begin{boxedminipage}{\linewidth}
\begin{alltt}\small

<?xml version="1.0" encoding="ISO-8859-1"?>
<stdn nodes="9" arcs="20" timeHorizon="44" name="stdn.xml">
    <node number="1">
        <penalty t="0" c1="1000" c2="1"/>
        <penalty t="1" c1="954" c2="47"/>
        ...
    </node>
    ...
    <node number="4">
        <wait t="0" time="1" c1="485" c2="629"/>
        <wait t="1" time="1" c1="491" c2="693"/>
        ...
    </node>
    ...
    <arc head="1" tail="2">
        <leavingTime t="0" c1="872" c2="22">
            <travelTime t="2" prob="250000"/>
            <travelTime t="3" prob="500000"/>
            <travelTime t="4" prob="250000"/>
        </leavingTime>
        ...
    </arc>
    ...
</stdn>
\end{alltt}
\end{boxedminipage}
  \caption{The xml output file.}
  \label{fig:xml}
\end{figure}

Output is per default written to $<$filename$>$.xml where $<$filename$>$ is the filename
specified by \textbf{-out}. The xml format is simple to understand and
illustrated in \Fref{fig:xml}. Note, the probabilities of the travel time
distribution is NOT normalized. This have to be done by the program which uses
the test instance.

In general, the xml format is very verbose resulting in large file sizes.
However, the xml file may be converted to a desired format, e.g.\ a
time-expanded hypergraph \citep{Pretolani98}, using an xslt stylesheet. For an
introduction to xml and xslt see \citet{Moller06}.
For very large test instances using a stylesheet may not be possible. Hence an option -f5
to generate the output in a time-expanded hypergraph format is also possible.

\endlatexonly

\htmlonly

\note{Change log}

\li 16-12-05: New version 1.5 finished with xml output.

\end{description}

\todo
   \li Check if functions and input can be made const.
   \li Check if input can be references.
   \li Check if some functions can be made private.

\endhtmlonly

*/

//------------------------------------------------------------------------------------

typedef int* IntPtr;



/** Class for representing a grid node and its forward star arcs.
\author Lars Relund Nielsen.
\version 1.66
*/
class GridNode
{
public:

    /** Constructor. Do not allocate memory which have to be done with AllocMem. */
    GridNode(){mem=false;}

    /** Deconstructor. Free memory automatically. */
    ~GridNode(){FreeMem();}

    /** Allocate memory for the gridNode.
    \param H Time-horizon.
    */
    void AllocMem(int H){c1 = new int[H+1][4];c2 = new int[H+1][4];mem=true;}

    /** Free memory. */
    void FreeMem(){if (mem) delete [] c1;delete [] c2;}

    int cWait1;     ///< First off-peak cost for waiting in the node.
    int cWait2;     ///< Second off-peak cost for waiting in the node.
    int traveltime[4]; ///< Mean off-peak travel times for grid arcs in forward star. Entrys are 0: north, 1: east, 2: south, 3: west.

    int (*c1)[4];   ///< First cost array for leaving at time t. Example c1[t][2] is the 1. cost for leaving the south arc at time t.
    int (*c2)[4];   ///< Second costs array leaving at time t. Example c2[t][2] is the 2. cost for leaving the south arc at time t.

private:
    bool mem;   ///< True if memory allocated.
};



//-----------------------------------------------------------------------------


/** The class for generating STD networks.

Assumes a underlying topological grid network of size \f$b\times h\f$. Every
arc is ``bi-directional'', but no arcs enter the origin and no arcs exit the
destination. Thus there are:

\li h(b-1)-2 arcs east.
\li b(h-1)-2 arcs south.
\li h(b-1) arcs west.
\li b(h-1) arcs north.

Grid nodes are numbered from the destination node 1 to the origin node
\f$b\cdot h\f$. This number increases by one south and by \f$h\f$ east. Thus,
if \f$u\f$ is the node with coordinates \f$(x,y)\f$ (\f$x\f$ right, \f$y\f$
down), then node \f$u\f$ have number \f[ (x-1)h + y \f] This number is used to
identify a node in the array of GridNode objects Moreover, each arc in a node
\f$u\f$ is numbered 0 (north), 1 (east), 2 (south) and 3 (west) if exists.

\author Lars Relund Nielsen.
\version 1.66
*/
class TegPeak
{
public:
    /** Run the generator and write the output to an xml file. */
    void Run();

    /** Read input parameters.
    Moreover, set the output file name and the flag for verbose, f5 output and xml output.
    */
    TegPeak(string filename,bool verbose,bool f5, bool xml);

    /** Free memory. */
    ~TegPeak();

private:

    /** Statistics class for internal statistics. */
    class Statistics {
    public:
        Statistics()
        {
            minT=minI=minC=INFIN;
            maxT=maxI=maxPdfSize=maxC=-1;
            avePdfSize=avePdfUb=0;
            ma=mh=hsize=0;
        }

        void Reset()
        {
            minT=minI=minC=INFIN;
            maxT=maxI=maxPdfSize=maxC=-1;
            avePdfSize=avePdfUb=0;
            ma=mh=hsize=0;
        }

        int minT,   ///< Minimum travel time used.
            maxT,   ///< Manimum travel time used.
            minI,   ///< Minimum interval travel time used.
            maxI,   ///< Maximum interval travel time used.
            maxPdfSize,     ///< Maximum size of a travel time distribution.
            maxC,   ///< Maximum travel cost used.
            minC,   ///< Minimum travel cost used.
            ma,     ///< Number of arcs in hypergraph.
            mh,     ///< Number of hyperarcs in hypergraph.
            hsize;  ///< Size of true hyperarcs (tails and heads).
        double avePdfSize,  ///< Average travel distribution size.
            avePdfUb;   ///< Average maximum travel time in all distributions.
    };

    /** Return the number of the node in the grid.
    \param x Base coordinate.
    \param y Height coordinate.
    */
    int GridNumber(int x,int y) {return (x-1)*itsHeight + y;}

    /** \brief Return the hypergraph node number.
     * \param gridN Grid number.
     * \param time Time instance.
     */
    int HNumber(int gridN,int time) {return (itsT-time-1)*itsDimension+gridN+1;}

    /** Given time t find the time it corresponds to in the first cycle. */
    int CycleTime(int t)
    {
        int cycle = t/itsCycleLength+1;
        t = t-(cycle-1)*itsCycleLength;  // time in first cycle
        return t;
    }

    /** Given reverse time t find the time instance. */
    int Time(int t)
    {
        return itsT-t+1;
    }

    /**  Check if a north arc exists for grid node gridNum. */
    bool North(int gridNum)
    {
        int x = 1 + (gridNum-1)/itsHeight;  // x koordinat
        int y = gridNum - itsHeight*(x-1);  // y koordinat
        if (y>1) return true;
        else return false;
    }

    /** Check if a east arc exists for grid node gridNum. */
    bool East(int gridNum)
    {
        int x = 1 + (gridNum-1)/itsHeight;  // x koordinat
        int y = gridNum - itsHeight*(x-1);  // y koordinat
        if (x<itsBase && gridNum!=itsDimension-itsHeight && gridNum!=1) return true;
        else return false;
    }

    /** Check if a south arc exists for grid node gridNum. */
    bool South(int gridNum)
    {
        int x = 1 + (gridNum-1)/itsHeight;  // x koordinat
        int y = gridNum - itsHeight*(x-1);  // y koordinat
        //cout << x << " " << y;
        if (y<itsHeight && gridNum!=itsDimension-1 && gridNum!=1) return true;
        else return false;
    }

    /** Check if a west arc exists for grid node gridNum. */
    bool West(int gridNum)
    {
        int x = 1 + (gridNum-1)/itsHeight;  // x koordinat
        int y = gridNum - itsHeight*(x-1);  // y koordinat
        if (x>1) return true;
        else return false;
    }

    /** Return the direction of the arc given arc (tail,head) in the grid.
    \return 0 (north), 1 (east), 2 (south), 3 (west).
    */
    int Direction(int head,int tail)
    {
        if (tail-1==head) return 0;
        if (tail+1==head) return 2;
        if (tail-itsHeight==head) return 3;
        if (tail+itsHeight==head) return 1;
        cout << "Error in Direction";
        exit(1);
    }

    /** Returns the head node of the grid arc with tail tail and direction dir. */
    int Head(int tail,int dir)
    {
        if (dir==0) return tail-1;
        if (dir==1) return tail+itsHeight;
        if (dir==2) return tail+1;
        if (dir==3) return tail-itsHeight;
        cout << "Error in Head";
        exit(1);
    }

    /** Round a double to an integer. */
    int Round(double i)
    {
        double j = i-(int)i;
        if (j>=0.5) return (int)i+1;
        else return (int)i;
    }

    /** Apply the random perturbation to a travel cost. */
    int Rand(int cost)
    {
        double c;
        if (itsRand==0) return cost;
        c = (double)cost*(1+(double)NumberGen.Int_number(-itsRand,itsRand)/1000);
        return MAX(1,Round(c));
    }

    /** Apply the random perturbation to a waiting cost. */
    int RandWait(int cost)
    {
        double c;
        if (itsRand==0) return cost;
        c = (double)cost*(1+(double)NumberGen.Int_number(-itsRand,itsRand)/1000);
        //cout << c << "\n";
        return MAX(1,Round(c));
    }


    /** Identify if we are in a peak.
    \return 0 if is in a non peak, 1 if is in the first transient part of the peak,
    2 if is in the pure peak part and 3 if is in the last transient part of the peak.
    */
    int Peak(int t);

    /** Generate mean travel times.
    Note to assure that the travel times for north and west arcs is the
    same with or without symmetry we sometimes generate costs not used.
    */
    void GenOffPeakMeanTraveltimes();

    /** Generate the waiting and travel costs.
    The costs are stored in the array itsGrid.
    Do not use the random element here applied later.
    \see GridNode.
    */
    void GenWeights();

    /** Generate waiting costs and times for nodes and leaving times.
    The generated costs and time is output to file.
    \note Waiting is always from t to t+1.
    */
    void GenWaiting();

    /** Generate the costs and travel time distributions.
    Write the generated costs and times to xml file if out is true.
    \pre GenOffPeakMeanTraveltimes and GenOffPeakWeights must have been called.
    \return The average max travel time for the travel time distributions.
    \param out True if output distributions and travel costs.
    */
    double GenTravel(bool out);

    /** Generate the cost and travel time distributions for arc (gTail,gHead).
    Write the generated cost and travel time distributions to xml file if out is true.
    \param mean Off-peak mean travel time.
    \param out True if output distributions and travel costs.
    \return The number travel time distributions generated.
    */
    int GenTravelArc(int gTail,int gHead,int mean,bool out);

    void GenerateGraf();

    /** Generate penalty costs for arrival to the destination and print to xml file.
    Note do always print the penalty costs to the xml file even if \f$ub_W=0\f$. In this
    case both penalty costs is just set to zero. This is done since then the xslt file for
    converting the output to a time-expanded hypergraph format can be done easiler.
    */
    void GenPenalties();       // i.e. generate arcs from `s' to (d,t) nodes

    /** Weight Pairs Generator for penalty costs for arrival to the destination at time t.
    \see GenPenalties.
    */
    void PenaltyW(int t);

    /** Generate travel costs according to the correlation type.
    The costs are stored in itsW1 and itsW2.
    */
    void WPair(int lb,int ub);

    /** Generate waiting costs according to the correlation type.
    The costs are stored in itsW1 and itsW2.
    */
    void WaitPair(int lb,int ub);


    void PrintFirstLineGraf();



    /** Print the input parameters. */
    void PrintInput();

    /** Check if input is okay. Exit program if not. */
    void CheckInput();

    /** Find the mean travel time at time t given the off-peak mean.
    */
    double FindMean(int t,int meanOffPeak);

    /** Find the cost at time t given the off-peak cost.
    We here assume that we increase the cost in peaks.*/
    double FindCostPlus(int t,int costOffPeak);

    /** Find the cost at time t given the off-peak cost.
    We here assume that we decrease the cost in peaks.*/
    double FindCostMinus(int t,int costOffPeak);

    /** Initialize the xml writer for the output. */
    void xmlInit();

    /** Free the xml writer and write the xml to a file. */
    void xmlFree();

    /** Start the element arc in xml file with attribute gHead and gTail.
    \param  gHead Head node in the grid.
    \param  gTail Tail node in the grid.
    */
    void startXmlArc(int gHead,int gTail);

    /** End the element arc in xml file. */
    void endXmlArc();

    /** Initialize the file for f5 output.
    Note that since the content of the first line depends on the next lines
    we first don't write the first line but afterwards insert it into the
    file (by creating a temporary file).
    */
    void f5Init();

    /** Close the f5 file. */
    void f5Free();

    /** Print error msg for xmlwriter */
    void PrintError(int rc) const
    {
        if(rc<0)
        {
            cout << "Error writing to file!\n"<<flush;
            exit(1);
        }
    }

    /** Initialize variables and random number generators. */
    void Initialize();

    /** Do preprocessing to find the average max travel time for each arc
    and leaving time.
    \return The average max travel time.
    */
    double Preprocess()
    {
        double ub;

        if (itsVerbose) cout << "\nDo some preprocessing to find time-horizon H ... ";
        ub = GenTravel(false);   // find an estimate of average max travel time in pdfs
        if (itsVerbose) cout << "done\n";
        return ub;
    }

    /** Allocate memory for the GridNodes.
    \pre itsT and itsDimension is set and memory not alread allocated.
    */
    void AllocMemGrid()
    {
        for(int i=0;i<itsDimension+1;i++) itsGrid[i].AllocMem(itsT);
    }

    /** Free memory for the GridNodes.
    \pre itsDimension is set.
    */
    void FreeMemGrid()
    {
        for(int i=0;i<itsDimension+1;i++) itsGrid[i].FreeMem();
    }

    /** Write the statistics to the xml file. */
    void WriteStat();

    /** Write the TEGP info to the xml file. */
    void WriteTegpInfo();

private:

    // input parameters
    int itsBase,             ///< Base of the grid
        itsHeight,           ///< Height of the grid.
        itsCycleLength,      ///< Time instances in a cycle.
        itsPeaks,            ///< Number of peak periodes in a cycle.
        itsPeakTrans,        ///< Time instances in the transient parts of the peak.
        itsPeakPure,         ///< Time instances in the pure part of the peak.
        itsPeakStart,        ///< Starting time for first peak.
        itsMeanIncrease,     ///< Increase in mean traveltime in peaks (pct)
        itsVarMean,          ///< Variance/mean ratio (pct).
        itsFSCmin,itsFSCmax, ///< Min and max cost in FS(s).
        itsTmin,itsTmax,     ///< Min and max cost on waiting arcs.
        itsTType,            ///< Travel time flag: 0 if all mean travel times are increased in peaks, 1 if north and west arcs do not have a peak effect, 2 if north and west arcs only have half the peak effect.
        itsStatic,           ///< Static flag: 0 if the network is not static in off-peaks, 1 if is static in off-peaks
        itsWaitMin,itsWaitMax,  ///< Min and max cost on waiting arcs.
        itsCmin,itsCmax,     ///< Min and max cost on harcs.
        itsSym,              ///< Symmetri flag. 1 if use symmetri, 0 otherwise.
        itsCorType,          ///< Correlation type.
        itsFS_s_Dep,         ///< FS(s) cost dependences.
        itsWaitDep,          ///< Waiting arcs cost dependences.
        itsHArcDep,          ///< Harcs cost dependences.
        itsRand,             ///< Random element (promille).
        itsSeed;

    // command line parameters
    string itsFilename;     ///< Filename for output.
    bool itsVerbose;        ///< Print verbose.
    bool itsF5;             ///< Save in f5 format.
    bool itsXml;            ///< Save in xml format.

    // variables for storing costs
    int itsW1,  ///< First cost.
        itsW2;  ///< Second cost.

    int itsDimension,           ///< Number of nodes in the grid (base*height).
        itsNonPeakLength,       ///< Length of off-peak
        itsFirstNonPeakLength,  ///< First off-peak length in the cycle.
        itsLastNonPeakLength,   ///< Last off-peak length in the cycle.
        itsPeakLength,          ///< Length of peak.
        itsT;                   ///< Number of timeinstances in total
    double itsPeriods;          ///< Number of cycles (can be fractional).

    IntPtr itsHArcArray,        ///< Contain tailnumbers.
        itsMultArray;           ///< Contain the multipliers.
    GridNode *itsGrid;          ///< Contain cost and traveltime info for each grid node.

    Random NumberGen;           ///< Random number routines.

    xmlTextWriterPtr pXmlWriter;    ///< Pointer to xml writer.
    ofstream f5File;                 ///< File for writing f5 format.

    Statistics itsStat;             ///< Internal statistics object.
};

//-----------------------------------------------------------------------------

#endif
