function _parse_neighborhood(tok::AbstractString)::Neighborhood
    if tok == "INC"
        return NB_INCLUSION
    elseif tok == "EXC"
        return NB_EXCLUSION
    elseif tok == "SYM_DIFF" || tok == "SYMDIFF"
        return NB_SYMDIFF
    else
        throw(ArgumentError("unknown neighborhood token: $tok"))
    end
end

function _neighborhood_token(nb::Neighborhood)::String
    nb == NB_INCLUSION && return "INC"
    nb == NB_EXCLUSION && return "EXC"
    return "SYM_DIFF"
end

function _parse_uncertainty(tok::AbstractString)::Uncertainty
    if tok == "NOMINAL"
        return U_NOMINAL
    elseif tok == "INTERVAL"
        return U_INTERVAL
    elseif tok == "CONT" || tok == "CONT_BUDGET"
        return U_CONT_BUDGET
    elseif tok == "DISC" || tok == "DISC_BUDGET"
        return U_DISC_BUDGET
    else
        throw(ArgumentError("unknown uncertainty token: $tok"))
    end
end

function _uncertainty_token(U::Uncertainty)::String
    U == U_NOMINAL && return "NOMINAL"
    U == U_INTERVAL && return "INTERVAL"
    U == U_CONT_BUDGET && return "CONT"
    return "DISC"
end

function _parse_params_line(line::AbstractString)::Params
    p = params_nominal()
    U = p.uncertainty
    nb = p.neighborhood
    k = p.k
    gamma = p.gamma
    delta = p.delta
    for tok in split(strip(line))
        kv = split(tok, '='; limit = 2)
        length(kv) == 2 || throw(ArgumentError("expected key=value in parameter line, got $tok"))
        key = kv[1]
        val = kv[2]
        if key == "neighborhood"
            nb = _parse_neighborhood(val)
        elseif key == "k"
            k = parse(Int32, val)
        elseif key == "uncertainty"
            U = _parse_uncertainty(val)
        elseif key == "gamma"
            gamma = parse(Float64, val)
        elseif key == "delta"
            delta = parse(Int32, val)
        else
            throw(ArgumentError("unknown parameter key: $key"))
        end
    end
    k < 0 && throw(ArgumentError("k < 0"))
    _check_budgets(gamma, delta)
    return Params(U, nb, k, gamma, delta)
end

function _parse_arc_line(line::AbstractString)
    toks = split(strip(line))
    length(toks) == 5 || throw(ArgumentError("truncated or malformed arc line: $line"))
    u = parse(Int32, toks[1])
    v = parse(Int32, toks[2])
    C = parse(Float64, toks[3])
    c_hat = parse(Float64, toks[4])
    d = parse(Float64, toks[5])
    return u, v, C, c_hat, d
end

"""
    parse_rrsp(filename) -> (Network, Params)

Read an RRSP instance from `filename`.

    s t
    neighborhood=INC k=0 uncertainty=NOMINAL gamma=0.0 delta=0
    u v C c_hat d
    ...

The parameter line is optional; if omitted, [`params_nominal`](@ref) is used.
Blank lines and comments (`# ...`) are skipped.
Vertex identifiers in the file become [`Graph`](@ref).`label`; dense indices
follow first-appearance order, starting with `s`, then `t` if distinct, then
arc endpoints in file order.

# Throws
- `ArgumentError`: if the file cannot be opened or does not conform to the
  format.
"""
function parse_rrsp(filename::AbstractString)::Tuple{Network,Params}
    isfile(filename) || throw(ArgumentError("cannot open RRSP instance file: $filename"))
    raw = readlines(filename)
    nlines = length(raw)
    idx = 1
    function peek_content()
        while idx <= nlines
            s = strip(raw[idx])
            if isempty(s) || startswith(s, "#")
                idx += 1
                continue
            end
            return s
        end
        return nothing
    end
    function take_content()
        s = peek_content()
        s === nothing && throw(ArgumentError("unexpected end of RRSP instance file: $filename"))
        idx += 1
        return s
    end

    peek_content() === nothing && throw(ArgumentError("empty RRSP instance file: $filename"))

    params = params_nominal()
    header = split(take_content())
    length(header) == 2 || throw(ArgumentError("header must be: s t"))
    s_uid = parse(Int32, header[1])
    t_uid = parse(Int32, header[2])
    nxt = peek_content()
    if nxt !== nothing && occursin('=', nxt)
        params = _parse_params_line(take_content())
    end
    arc_rows = Tuple{Int32,Int32,Float64,Float64,Float64}[]
    while true
        row = peek_content()
        row === nothing && break
        push!(arc_rows, _parse_arc_line(take_content()))
    end

    labels_order = Int32[]
    uid_to_dense = Dict{Int32,Int32}()
    function dense_of(uid::Int32)::Int32
        if !haskey(uid_to_dense, uid)
            push!(labels_order, uid)
            uid_to_dense[uid] = Int32(length(labels_order))
        end
        return uid_to_dense[uid]
    end
    dense_of(s_uid)
    dense_of(t_uid)
    m = length(arc_rows)
    tail = Vector{Int32}(undef, m)
    head = Vector{Int32}(undef, m)
    C = Vector{Float64}(undef, m)
    c_hat = Vector{Float64}(undef, m)
    d = Vector{Float64}(undef, m)
    for a in 1:m
        u, v, Ca, cha, da = arc_rows[a]
        tail[a] = dense_of(u)
        head[a] = dense_of(v)
        C[a] = Ca
        c_hat[a] = cha
        d[a] = da
    end
    n = length(labels_order)
    g = graph_new(n, tail, head; label = labels_order)
    costs = costs_new(C, c_hat, d)
    net = network_new(g, uid_to_dense[s_uid], uid_to_dense[t_uid], costs)
    return net, params
end

"""
    write_rrsp(filename, net, params)

Write `net` and `params` in the file format of [`parse_rrsp`](@ref).
"""
function write_rrsp(filename::AbstractString, net::Network, params::Params)
    open(filename, "w") do io
        g = net.graph
        println(io, "$(g.label[net.s]) $(g.label[net.t])")
        nb = _neighborhood_token(params.neighborhood)
        U = _uncertainty_token(params.uncertainty)
        println(
            io,
            "neighborhood=$nb k=$(Int(params.k)) uncertainty=$U gamma=$(params.gamma) delta=$(Int(params.delta))",
        )
        for a in 1:Int(g.m)
            println(
                io,
                "$(g.label[g.tail[a]]) $(g.label[g.head[a]]) $(net.costs.C[a]) $(net.costs.c_hat[a]) $(net.costs.d[a])",
            )
        end
    end
    return nothing
end

write_rrsp(filename::AbstractString, net::Network) = write_rrsp(filename, net, params_nominal())
