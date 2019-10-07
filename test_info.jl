using Coverage, OrdinaryDiffEq
info = analyze_malloc(dirname(pathof(OrdinaryDiffEq)))
info = info[findall(x->x.bytes != 0, info)]
