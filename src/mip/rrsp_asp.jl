function _asp_parents(tree::AspTree)::Tuple{Vector{Int},Vector{UInt8}}
    nn = Int(tree.n_nodes)
    parent = zeros(Int, nn)
    isleft = zeros(UInt8, nn)
    @inbounds for i in 1:nn
        nd = tree.nodes[i]
        nd.op == ASP_LEAF && continue
        lc = Int(nd.left)
        rc = Int(nd.right)
        parent[lc] = i
        parent[rc] = i
        isleft[lc] = 0x01
    end
    return parent, isleft
end

function _asp_helper(
    i::Int,
    l::Int,
    root::Int,
    parent::Vector{Int},
    isleft::Vector{UInt8},
    nodes::Vector{AspNode},
    DualRoot,
    DualParL,
    DualParR,
    DualSer,
    L::Int,
)
    if i == root
        return DualRoot[l]
    end
    p = parent[i]
    if nodes[p].op == ASP_PARALLEL
        return isleft[i] == 0x01 ? DualParL[p, l] : DualParR[p, l]
    end
    expr = JuMP.AffExpr(0.0)
    if isleft[i] == 0x01
        @inbounds for j in 0:(L - l)
            JuMP.add_to_expression!(expr, DualSer[p, l + j, j])
        end
    else
        @inbounds for tot in l:L
            JuMP.add_to_expression!(expr, DualSer[p, tot, l])
        end
    end
    return expr
end

function _solve_rrsp_cont_asp(
    net::Network,
    tree::AspTree,
    params::Params,
    solver::Solver,
    t0::UInt64,
)::Solution
    params.neighborhood == NB_INCLUSION || return _solve_rrsp_cont_bg(net, params, solver, t0)
    g = net.graph
    m = Int(g.m)
    nn = Int(tree.n_nodes)
    L = Int(params.k)
    C = net.costs.C
    c_hat = net.costs.c_hat
    ddev = net.costs.d
    gamma = params.gamma
    root = Int(tree.root)
    parent, isleft = _asp_parents(tree)
    nodes = tree.nodes
    bigM = _sum_vec(c_hat) + _sum_vec(ddev) + 1.0
    model = _new_mip(solver)
    model === nothing && return _not_impl(m, METHOD_MIP, t0)
    xnode = JuMP.@variable(model, [1:nn], Bin)
    DualRoot = JuMP.@variable(model, [0:L], lower_bound = 0.0)
    DualLeaf = JuMP.@variable(model, [1:m, 0:L], lower_bound = 0.0)
    DualParL = JuMP.@variable(model, [1:nn, 0:L], lower_bound = 0.0)
    DualParR = JuMP.@variable(model, [1:nn, 0:L], lower_bound = 0.0)
    DualSer = JuMP.@variable(model, [1:nn, 0:L, 0:L], lower_bound = 0.0)
    beta = JuMP.@variable(model, [1:m], lower_bound = 0.0)
    theta = JuMP.@variable(model, lower_bound = 0.0)
    w0 = JuMP.@variable(model, [1:m], lower_bound = 0.0)
    w1 = JuMP.@variable(model, [1:m], lower_bound = 0.0)
    JuMP.@constraint(model, xnode[root] == 1)
    @inbounds for i in 1:nn
        nd = nodes[i]
        nd.op == ASP_LEAF && continue
        lc = Int(nd.left)
        rc = Int(nd.right)
        if nd.op == ASP_SERIES
            JuMP.@constraint(model, xnode[lc] == xnode[i])
            JuMP.@constraint(model, xnode[rc] == xnode[i])
        else
            JuMP.@constraint(model, xnode[lc] + xnode[rc] == xnode[i])
        end
    end
    JuMP.@constraint(model, sum(DualRoot[l] for l in 0:L) >= 1)
    @inbounds for i in 1:nn
        nd = nodes[i]
        for l in 0:L
            h = _asp_helper(i, l, root, parent, isleft, nodes, DualRoot, DualParL, DualParR, DualSer, L)
            if nd.op == ASP_LEAF
                a = Int(nd.arc)
                JuMP.@constraint(model, DualLeaf[a, l] >= h)
            elseif nd.op == ASP_SERIES
                JuMP.@constraint(model, sum(DualSer[i, l, j] for j in 0:l) >= h)
            else
                JuMP.@constraint(model, DualParL[i, l] + DualParR[i, l] >= h)
            end
        end
    end
    @inbounds for a in 1:m
        if L >= 1
            JuMP.@constraint(model, beta[a] + theta >= DualLeaf[a, 0] + DualLeaf[a, 1])
        else
            JuMP.@constraint(model, beta[a] + theta >= DualLeaf[a, 0])
        end
        JuMP.@constraint(model, w0[a] <= xnode[a])
        JuMP.@constraint(model, w0[a] <= DualLeaf[a, 0])
        JuMP.@constraint(model, w0[a] >= DualLeaf[a, 0] - (1.0 - xnode[a]))
        if L >= 1
            JuMP.@constraint(model, w1[a] <= xnode[a])
            JuMP.@constraint(model, w1[a] <= DualLeaf[a, 1])
            JuMP.@constraint(model, w1[a] >= DualLeaf[a, 1] - (1.0 - xnode[a]))
        end
    end
    dummy = JuMP.AffExpr(0.0)
    if L >= 2
        @inbounds for a in 1:m
            for l in 2:L
                JuMP.add_to_expression!(dummy, DualLeaf[a, l])
            end
        end
    end
    linM = JuMP.AffExpr(0.0)
    @inbounds for a in 1:m
        JuMP.add_to_expression!(linM, DualLeaf[a, 0])
        JuMP.add_to_expression!(linM, -1.0, w0[a])
        if L >= 1
            JuMP.add_to_expression!(linM, w1[a])
        end
    end
    c_hat_term = JuMP.AffExpr(0.0)
    @inbounds for a in 1:m
        JuMP.add_to_expression!(c_hat_term, c_hat[a], DualLeaf[a, 0])
        if L >= 1
            JuMP.add_to_expression!(c_hat_term, c_hat[a], DualLeaf[a, 1])
        end
    end
    JuMP.@objective(
        model,
        Min,
        sum(C[a] * xnode[a] for a in 1:m) + c_hat_term + bigM * linM + bigM * dummy +
        sum(ddev[a] * beta[a] for a in 1:m) + gamma * theta,
    )
    st, gap = _optimize_mip(model)
    dt = (time_ns() - t0) / 1e9
    if !_mip_usable(model, st)
        return _with_mip_diag(
            solution_empty(m; status = st, method = METHOD_MIP, time_sec = dt, mip_gap = gap),
            model,
        )
    end
    chi = _chi_from_vars(xnode, m)
    xp = try
        path_from_chi(g, net.s, net.t, chi)
    catch e
        e isa ArgumentError || rethrow()
        return _with_mip_diag(
            solution_empty(m; status = ST_ERROR, method = METHOD_MIP, time_sec = dt),
            model,
        )
    end
    z_mip = Float64(JuMP.objective_value(model))
    return _with_mip_diag(_rrsp_from_x(net, params, xp, z_mip, gap, st, solver, t0), model)
end
