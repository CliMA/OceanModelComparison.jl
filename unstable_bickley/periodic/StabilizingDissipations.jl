module StabilizingDissipations

using ClimateMachine.Mesh.Grids: min_node_distance

struct StabilizingDissipation{T}
    κh_min :: T 
    κh_max :: T
    νh_min :: T
    νh_max :: T
    smoothness_exponent :: T
    minimum_node_distance :: T
    Δu :: T
    Δθ :: T
end

function effective_grid_spacing(grid)
    return 1
end

function StabilizingDissipation(grid; 
                                minimum_node_distance = min_node_distance(grid, HorizontalDirection()),
                                Δu = 1,
                                Δθ = 1,
                                κh_max = 0,
                                κh_min = 0,
                                νh_max = 0,
                                νh_min = 0,
                                smoothness_exponent = 2)

    FT = eltype(grid)

    return StabilizingDissipation(
                                  FT(κh_min),
                                  FT(κh_max),
                                  FT(νh_min),
                                  FT(νh_max),
                                  FT(smoothness_exponent),
                                  FT(minimum_node_distance),
end


end # module
