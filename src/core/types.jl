# Vertices are identified with the integers 1, …, n.
# Arcs are identified with the integers 1, …, m.

"""
    Neighborhood

Neighborhood used to restrict recovery from a committed first-stage path.
The integer parameter ``k`` of [`Params`](@ref) is the neighborhood size.

# Values
- `NB_INCLUSION`: admissible recoveries ``y`` satisfy ``|y \\setminus x| \\leq k``.
- `NB_EXCLUSION`: admissible recoveries ``y`` satisfy ``|x \\setminus y| \\leq k``.
- `NB_SYMDIFF`: admissible recoveries ``y`` satisfy ``|x \\triangle y| \\leq k``.
"""
@enum Neighborhood::Int32 begin
    NB_INCLUSION = 0
    NB_EXCLUSION = 1
    NB_SYMDIFF   = 2
end

"""
    Uncertainty

Second-stage cost uncertainty set.

# Values
- `U_NOMINAL`: the singleton ``\\{\\hat{c}\\}``. The budget fields of [`Params`](@ref)
  are ignored.
- `U_INTERVAL`: each arc ``a`` may take any cost in ``[\\hat{c}_a, \\hat{c}_a + d_a]``.
  The budget fields of [`Params`](@ref) are ignored.
- `U_CONT_BUDGET`: interval uncertainty with a continuous budget ``Γ``
  ([`Params`](@ref).`gamma`).
- `U_DISC_BUDGET`: interval uncertainty with a discrete budget ``Δ``
  ([`Params`](@ref).`delta`).
"""
@enum Uncertainty::Int32 begin
    U_NOMINAL     = 0
    U_INTERVAL    = 1
    U_CONT_BUDGET = 2
    U_DISC_BUDGET = 3
end

"""
    Method

Algorithmic method requested of a solver.

# Values
- `METHOD_AUTO`: use a combinatorial algorithm when the case is known to be
  polynomial-time solvable, and a mixed-integer formulation otherwise.
- `METHOD_COMB`: combinatorial algorithm only; return `ST_NOT_IMPL` if none
  is available.
- `METHOD_MIP`: mixed-integer (or linear) formulation only.
"""
@enum Method::Int32 begin
    METHOD_AUTO = 0
    METHOD_COMB = 1
    METHOD_MIP  = 2
end

"""
    Status

Completion status of a solver or evaluation call.

# Values
- `ST_OK`: an optimal solution was obtained.
- `ST_INFEASIBLE`: no feasible solution exists (in particular, no ``s``–``t`` path).
- `ST_TIME_LIMIT`: the solver reached its time limit.
- `ST_UNBOUNDED`: the formulation was reported unbounded (should not occur
  under nonnegative costs).
- `ST_NOT_IMPL`: the requested method is not implemented for this case.
- `ST_ERROR`: the call failed for a reason not covered above.
"""
@enum Status::Int32 begin
    ST_OK         = 0
    ST_INFEASIBLE = 1
    ST_TIME_LIMIT = 2
    ST_UNBOUNDED  = 3
    ST_NOT_IMPL   = 4
    ST_ERROR      = 5
end

"""
    AspOp

Composition operation at a node of an arc-series-parallel decomposition tree.

# Values
- `ASP_LEAF`: the node represents a single arc.
- `ASP_SERIES`: series composition of the two children.
- `ASP_PARALLEL`: parallel composition of the two children.
"""
@enum AspOp::Int32 begin
    ASP_LEAF     = 0
    ASP_SERIES   = 1
    ASP_PARALLEL = 2
end

"""
    Graph

A directed multigraph ``G = (V, A)`` in forward- and reverse-star form.

Vertices are the integers `1, …, n` and arcs the integers `1, …, m`. An arc is
identified by its index, not by the pair of its endpoints; parallel arcs are
therefore permitted.

The arrays `tail` and `head` store the endpoints of each arc: arc `a` leaves
`tail[a]` and enters `head[a]`.

The forward star of vertex `v` is the block of arc indices
`out_arc[first_out[v]:first_out[v+1]-1]`, i.e. all arcs with tail `v`.
The reverse star is `in_arc[first_in[v]:first_in[v+1]-1]`, all arcs with
head `v`. The arrays `first_out` and `first_in` have length `n+1` and satisfy
`first_out[1] = 1`, `first_out[n+1] = m+1`, and likewise for `first_in`.
They are offsets into `out_arc` and `in_arc`, not lists of vertices.

The vector `label` stores input identifiers for the vertices. Algorithms
address vertices exclusively by the dense indices `1, …, n`.

Graphs are constructed by [`graph_new`](@ref) and are not mutated thereafter.

# Fields
- `n`: number of vertices, ``|V|``.
- `m`: number of arcs, ``|A|``.
- `tail`: `tail[a]` is the tail of arc `a`.
- `head`: `head[a]` is the head of arc `a`.
- `first_out`: `first_out[v]` is the first index in `out_arc` of the forward star of `v`.
- `out_arc`: arc indices grouped by tail, in increasing order of arc index.
- `first_in`: `first_in[v]` is the first index in `in_arc` of the reverse star of `v`.
- `in_arc`: arc indices grouped by head, in increasing order of arc index.
- `label`: external identifier of each vertex.
"""
struct Graph
    n::Int32
    m::Int32
    tail::Vector{Int32}
    head::Vector{Int32}
    first_out::Vector{Int32}
    out_arc::Vector{Int32}
    first_in::Vector{Int32}
    in_arc::Vector{Int32}
    label::Vector{Int32}
end

"""
    Costs

Arc costs of a two-stage instance, stored as arrays parallel to the arc index
of a [`Graph`](@ref).

# Fields
- `C`: first-stage cost ``C_a`` of each arc.
- `c_hat`: nominal second-stage cost ``\\hat{c}_a`` of each arc.
- `d`: maximum second-stage deviation ``d_a \\geq 0`` of each arc.
"""
struct Costs
    C::Vector{Float64}
    c_hat::Vector{Float64}
    d::Vector{Float64}
end

"""
    Network

An immutable network instance: a digraph, a pair of terminals, and arc costs.

A network does not include a neighborhood, a recovery parameter, or an
uncertainty budget. Those quantities belong to [`Params`](@ref).

# Fields
- `graph`: the underlying multidigraph.
- `s`: origin vertex.
- `t`: destination vertex.
- `costs`: first-stage and second-stage arc costs.
"""
struct Network
    graph::Graph
    s::Int32
    t::Int32
    costs::Costs
end

"""
    Params

Parameters of a robust or recoverable-robust problem on a fixed [`Network`](@ref).

This is the object varied in experimental sweeps. The field `gamma` is the
continuous budget ``Γ`` and is used only when `uncertainty == U_CONT_BUDGET`.
The field `delta` is the discrete budget ``Δ`` and is used only when
`uncertainty == U_DISC_BUDGET`.

# Fields
- `uncertainty`: second-stage uncertainty set.
- `neighborhood`: recovery neighborhood.
- `k`: neighborhood size ``k``.
- `gamma`: continuous budget ``Γ``.
- `delta`: discrete budget ``Δ``.
"""
struct Params
    uncertainty::Uncertainty
    neighborhood::Neighborhood
    k::Int32
    gamma::Float64
    delta::Int32
end

"""
    Path

A (simple) ``s``–``t`` path in a [`Graph`](@ref), stored in two equivalent forms.

The characteristic vector `chi` has length `m` and satisfies `chi[a] = 1` if and
only if arc `a` lies on the path. The vector `seq` lists the same arcs in
traversal order from `s` to `t`. After construction by [`path_from_seq`](@ref)
or [`path_from_chi`](@ref), the two representations are required to be consistent.

An empty sequence together with a zero characteristic vector denotes the
absence of a path (no feasible recovery, unused stage, or infeasible instance),
or the trivial path when ``s = t``.

# Fields
- `chi`: characteristic vector of the path on the arc set, entries in ``\\{0,1\\}``.
- `seq`: arc indices in the order of traversal; empty if no path is represented.
"""
struct Path
    chi::Vector{UInt8}
    seq::Vector{Int32}
end

"""
    Solution

Result of a solver or evaluation call.

The scalar `z` is the objective of the problem that was solved (nominal, robust,
incremental, recoverable, adversarial, or recoverable-robust). Nominal and
sampled evaluations of a committed path are obtained from the evaluation API,
not by overloading `z`.

If `status != ST_OK`, the path fields are empty and the cost fields are `Inf`.

# Fields
- `first`: first-stage path (or the unique path, for single-stage problems).
- `second`: second-stage (recovery) path; empty when the call does not produce one.
- `z`: objective value of the call.
- `z_first`: first-stage cost ``C(x)`` of `first`.
- `z_second`: second-stage cost of `second` under the scenario used by the call.
- `status`: completion status.
- `method_used`: method that produced the result (`METHOD_COMB` or `METHOD_MIP`).
- `time_sec`: wall-clock time of the call, in seconds.
- `mip_gap`: relative MIP gap; `0` when the method is combinatorial.
- `mip_nodes`: branch-and-bound nodes; `0` when combinatorial or unknown.
- `n_binaries`: binary variables in the formulation; `0` when combinatorial.
- `n_constraints`: JuMP constraints counted after build; `0` when combinatorial.
"""
struct Solution
    first::Path
    second::Path
    z::Float64
    z_first::Float64
    z_second::Float64
    status::Status
    method_used::Method
    time_sec::Float64
    mip_gap::Float64
    mip_nodes::Int32
    n_binaries::Int32
    n_constraints::Int32
end

function Solution(
    first::Path,
    second::Path,
    z::Float64,
    z_first::Float64,
    z_second::Float64,
    status::Status,
    method_used::Method,
    time_sec::Float64,
    mip_gap::Float64,
)
    return Solution(
        first,
        second,
        z,
        z_first,
        z_second,
        status,
        method_used,
        time_sec,
        mip_gap,
        Int32(0),
        Int32(0),
        Int32(0),
    )
end

"""
    Solver

Configuration of a solution method, passed to every `solve_*` call.

There is no process-wide optimizer. The field `optimizer` is a JuMP-compatible
optimizer constructor, or `nothing` when only combinatorial methods are to be
used.

# Fields
- `optimizer`: JuMP optimizer constructor, or `nothing`.
- `method`: requested method ([`Rrsp.Method`](@ref)).
- `time_limit`: time limit in seconds; `Inf` for none.
- `mip_gap`: relative MIP optimality tolerance.
- `threads`: thread limit for the MIP solver; `0` leaves it to the solver.
- `silent`: if `true`, suppress solver output.
"""
struct Solver
    optimizer::Any
    method::Method
    time_limit::Float64
    mip_gap::Float64
    threads::Int32
    silent::Bool
end

"""
    AspNode

A node of an arc-series-parallel decomposition tree.

A leaf stores the arc it represents in `arc` and has `left = right = 0`. An
internal node stores the indices of its children in `left` and `right` and has
`arc = 0`. The fields `s` and `t` are the terminals of the subgraph composed at
this node.

# Fields
- `op`: composition operation ([`AspOp`](@ref)).
- `left`: index of the left child, or `0`.
- `right`: index of the right child, or `0`.
- `arc`: arc index if `op == ASP_LEAF`, otherwise `0`.
- `s`: origin terminal of the subgraph.
- `t`: destination terminal of the subgraph.
"""
struct AspNode
    op::AspOp
    left::Int32
    right::Int32
    arc::Int32
    s::Int32
    t::Int32
end

"""
    AspTree

Binary decomposition tree of an arc-series-parallel digraph, as returned by
[`asp_decompose`](@ref).

Nodes are addressed by the integers `1, …, n_nodes`. Child indices in
[`AspNode`](@ref) refer to this numbering.

# Fields
- `nodes`: nodes of the tree.
- `root`: index of the root.
- `n_nodes`: number of nodes, equal to `2m-1` for a graph with `m` arcs.
"""
struct AspTree
    nodes::Vector{AspNode}
    root::Int32
    n_nodes::Int32
end
