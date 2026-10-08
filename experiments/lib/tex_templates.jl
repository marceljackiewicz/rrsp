# Fixed text of the LaTeX tables, copied from the thesis chapter so that the
# regenerated tables are identical in everything except the numbers.

const TABLE_INTERVAL_INSTANCES_HEAD = raw"""
\begin{table}[h!]
    \centering
    \caption{
        Interval instances.
        Here $\ell$ is the minimum number of arcs of an $\STPath$ path, and the last column is the maximum.
        The random-DAG row is one draw.
    }
    \label{tab:experiments-interval-instances}
    \begin{tabular}{lrrrr}
        \hline
        Family & $|\VerticesSet|$ & $|\ArcsSet|$ & $\ell$ & maximum number of arcs \\
        \hline
"""
const TABLE_INTERVAL_INSTANCES_TAIL = raw"""
        \hline
    \end{tabular}
\end{table}
"""

const TABLE_INTERVAL_VOR_HEAD = raw"""
\begin{table}[h!]
\centering
\caption{
    Mean value of recovery (\%) shown in Figure~\ref{fig:scale-vor-independent}.
    Column \(j\) corresponds to the horizontal coordinate \(j/12\), for $0 \leq j \leq 12$.
    For the ASP row that coordinate is \(\NeighborhoodSize/(3\ell)\).
    For the random DAG, intermediate columns are the linearly interpolated means.
}
\label{tab:experiments-interval-vor}
\footnotesize
\setlength{\tabcolsep}{3.5pt}
\begin{tabular}{l*{13}{r}}
\hline
\(12\cdot\NeighborhoodSize/\ell\) & $0$ & $1$ & $2$ & $3$ & $4$ & $5$ & $6$ & $7$ & $8$ & $9$ & $10$ & $11$ & $12$ \\
\hline
"""
const TABLE_INTERVAL_VOR_TAIL = raw"""
\hline
\end{tabular}
\end{table}
"""

const TABLE_BOTTLENECK_VOR_HEAD = raw"""
\begin{table}[h!]
\centering
\caption{
    Value of recovery (\%) shown in Figure~\ref{fig:scale-vor-bottleneck}.
    Column \(j\) reports the value at horizontal coordinate \(j/12\).
    For the layered digraph, the grid, and the random DAG this is \(\NeighborhoodSize/\ell = j/12\).
    For the ASP row it is \(\NeighborhoodSize/(3\ell) = j/12\).
    For the random DAG, columns \(3\), \(7\), and \(11\) are obtained by linear interpolation.
}
\label{tab:experiments-bottleneck-vor}
\footnotesize
\setlength{\tabcolsep}{3.5pt}
\begin{tabular}{l*{13}{r}}
\hline
\(12\cdot\NeighborhoodSize/\ell\) & $0$ & $1$ & $2$ & $3$ & $4$ & $5$ & $6$ & $7$ & $8$ & $9$ & $10$ & $11$ & $12$ \\
\hline
"""
const TABLE_BOTTLENECK_VOR_TAIL = raw"""
\hline
\end{tabular}
\end{table}
"""

const TABLE_CONT_INSTANCES_HEAD = raw"""
\begin{table}[h!]
    \centering
    \caption{
        Continuous-budget instances.
        Here $\ell$ is the minimum number of arcs of an $\STPath$ path, and the last column is the maximum.
        The random-DAG row is the range over the ten draws.
    }
    \label{tab:experiments-cont-instances}
    \begin{tabular}{lrrrrr}
        \hline
        Family & $|\VerticesSet|$ & $|\ArcsSet|$ & number of $\STPath$ paths & $\ell$ & maximum number of arcs \\
        \hline
"""
const TABLE_CONT_INSTANCES_TAIL = raw"""
        \hline
    \end{tabular}
\end{table}
"""

const TABLE_CONT_VOR_HEAD = raw"""
\begin{table}[h!]
\centering
\caption{
    Mean value of recovery (\%) shown in Figure~\ref{fig:scale-vor-cont}.
    Column \(j\) corresponds to the horizontal coordinate \(j/12\), for $0 \leq j \leq 12$.
    For the ASP rows that coordinate is \(\NeighborhoodSize/(3\ell)\).
    For the random DAG, intermediate columns are the linearly interpolated means.
    Fractions $0$ and $1$ are identically zero and are omitted.
}
\label{tab:experiments-cont-vor}
\scriptsize
\setlength{\tabcolsep}{2.6pt}
\begin{tabular}{cl*{13}{r}}
\hline
Fraction & Family & $0$ & $1$ & $2$ & $3$ & $4$ & $5$ & $6$ & $7$ & $8$ & $9$ & $10$ & $11$ & $12$ \\
\hline
"""
const TABLE_CONT_VOR_TAIL = raw"""
\hline
\end{tabular}
\end{table}
"""

const TABLE_DISC_INSTANCES_HEAD = raw"""
\begin{table}[h!]
    \centering
    \caption{
        Discrete-budget instances.
        Here $\ell$ is the minimum number of arcs of an $\STPath$ path, and the last column is the maximum.
        The random-DAG row is the range over the ten draws.
    }
    \label{tab:experiments-disc-instances}
    \begin{tabular}{lrrrrr}
        \hline
        Family & $|\VerticesSet|$ & $|\ArcsSet|$ & number of $\STPath$ paths & $\ell$ & maximum number of arcs \\
        \hline
"""
const TABLE_DISC_INSTANCES_TAIL = raw"""
        \hline
    \end{tabular}
\end{table}
"""

const TABLE_DISC_VOR_HEAD = raw"""
\begin{table}[h!]
\centering
\caption{
    Mean value of recovery (\%) shown in Figure~\ref{fig:scale-vor-disc}.
    The columns are analogous to those of Table~\ref{tab:experiments-cont-vor}.
    On the layered digraph and the grid, $\ell = 8$, so a column falling between solved neighborhood sizes is linearly interpolated.
}
\label{tab:experiments-disc-vor}
\scriptsize
\setlength{\tabcolsep}{2.6pt}
\begin{tabular}{cl*{13}{r}}
\hline
Fraction & Family & $0$ & $1$ & $2$ & $3$ & $4$ & $5$ & $6$ & $7$ & $8$ & $9$ & $10$ & $11$ & $12$ \\
\hline
"""
const TABLE_DISC_VOR_TAIL = raw"""
\hline
\end{tabular}
\end{table}
"""

# Figure captions and labels (end of each figure environment), also copied from the thesis.
const FIGURE_INTERVAL_TAIL = raw"""
\caption{Mean value of recovery under independent costs,
$\FirstStageCostVector = \NominalSecondStageCostVector$,
and the inclusion neighborhood.
The layered, grid, and arc-series-parallel curves each average $30$ cost draws on one digraph.
The random-DAG curve averages ten digraphs, with three cost draws on each.
The arc-series-parallel series is drawn against $\NeighborhoodSize/(3\ell)$.}
\label{fig:scale-vor-independent}
\end{center}
\end{figure}
"""

const FIGURE_BOTTLENECK_TAIL = raw"""
\caption{Value of recovery for the designated-path costs,
with $\FirstStageCostVector = \NominalSecondStageCostVector$
and the inclusion neighborhood.
One $\STPath$ has nominal cost $1$ and deviation $100$ on each arc;
every other arc has nominal cost $8$ and deviation $1$.
The arc-series-parallel series is drawn against $\NeighborhoodSize/(3\ell)$.}
\label{fig:scale-vor-bottleneck}
\end{center}
\end{figure}
"""

const FIGURE_CONT_TAIL = raw"""
\caption{Mean value of recovery under continuous budgeted uncertainty,
$\FirstStageCostVector = 0$, and the inclusion neighborhood.
Each curve is a fraction of the least budget at which some scenario gives every $\STPath$ path cost at least equal to the optimal value of $\RobProblemShort$ under $\IntervalUncertaintySet$.
Fractions $0$ and $1$ are identically zero and are omitted.
The layered, grid, and arc-series-parallel curves each average $30$ cost draws on one digraph.
The random-DAG curve averages ten digraphs, with three cost draws on each.
The arc-series-parallel series is drawn against $\NeighborhoodSize/(3\ell)$.}
\label{fig:scale-vor-cont}
\end{center}
\end{figure}
"""

const FIGURE_DISC_TAIL = raw"""
\caption{Mean value of recovery under discrete-budgeted uncertainty,
$\FirstStageCostVector = 0$, and the inclusion neighborhood,
averaged as in Figure~\ref{fig:scale-vor-cont}.
A fraction lying between two integer budgets is the linear interpolation of the optimal values.}
\label{fig:scale-vor-disc}
\end{center}
\end{figure}
"""

