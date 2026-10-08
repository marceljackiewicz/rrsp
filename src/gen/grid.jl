"""
    gen_grid(nrows, ncols) -> Graph

Directed grid DAG: `nrows × ncols` vertices, arcs only rightward and
downward. Origin is vertex `1` (top-left); destination is `nrows*ncols`
(bottom-right). Every ``s``–``t`` path has
`(nrows-1)+(ncols-1)` arcs.
"""
function gen_grid(nrows::Integer, ncols::Integer)::Graph
    nrows >= 1 || throw(ArgumentError("nrows < 1"))
    ncols >= 1 || throw(ArgumentError("ncols < 1"))
    nr = Int(nrows)
    nc = Int(ncols)
    tail = Int32[]
    head = Int32[]
    @inline vid(r, c) = Int32((r - 1) * nc + c)
    for r in 1:nr
        for c in 1:nc
            u = vid(r, c)
            if c < nc
                push!(tail, u)
                push!(head, vid(r, c + 1))
            end
            if r < nr
                push!(tail, u)
                push!(head, vid(r + 1, c))
            end
        end
    end
    return graph_new(nr * nc, tail, head)
end
