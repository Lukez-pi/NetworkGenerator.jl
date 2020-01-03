module NetworkGenerator
using StatsBase

RXN_MECH = ["UNIUNI", "UNIBI", "BIUNI", "BIBI"]
RXN_MECH_WEIGHT = [0.1, 0.2, 0.3, 0.4]

export randomize_rxns
export generate_model
export find_boundary_species
export print_dict

function print_dict(d::Dict, pre=1)
    for (k,v) in d
        if typeof(v) <: Dict
            s = "$(repr(k)) => "
            println(join(fill(" ", pre)) * s)
            print_dict(v, pre+1+length(s))
        else
            println(join(fill(" ", pre)) * "$(repr(k)) => $(repr(v))")
        end
    end
    nothing
end

function pick_reaction()
    if sum(RXN_MECH_WEIGHT) != 1
        error("the sum of weights don't add up to 1")
    end
    return sample(RXN_MECH, Weights(RXN_MECH_WEIGHT))
end

function randomize_rxns(species::Vector{String}, geneSpecies::Vector{String}, nRxns::Int64; rxn_specs=Dict(), reacting_species = Set{String}(), delete_rxns = Dict())
    addition_specs = Dict()
    reacting_species_init = copy(reacting_species)
    missing_species = true
    derivedSpecies = setdiff(Set(species), Set(geneSpecies))
    while missing_species
        nGeneratedRxns = 0
        while nGeneratedRxns < nRxns
            conflicting_rxns = true
            while conflicting_rxns
                rxn_mech = pick_reaction()
                if rxn_mech == "UNIUNI"
                    rct, prd, cat = sample(collect(species), 3, replace=false)
                    rct = [rct]
                    prd = [prd]
                    cat = [cat]
                elseif rxn_mech == "UNIBI"
                    rct, prd1, prd2 = sample(collect(species), 3, replace=false)
                    rct = [rct]
                    prd = [prd1, prd2]
                    cat = []
                elseif rxn_mech == "BIUNI"
                    rct1, rct2, prd = sample(collect(species), 3, replace=false)
                    rct = [rct1, rct2]
                    prd = [prd]
                    cat = []
                elseif rxn_mech == "BIBI"
                    rct1, rct2, prd1, prd2 = sample(collect(species), 4, replace=false)
                    rct = [rct1, rct2]
                    prd = [prd1, prd2]
                    cat = []
                else
                    error("Invalid reaction mechanism: ", rxn_mech)
                end

                forward_key, reverse_key = convert_rctprd_to_keys(rct, prd)
                if haskey(addition_specs, reverse_key)
                    if addition_specs[reverse_key]["cat"] == cat
                        continue
                    end
                end
                if haskey(rxn_specs, reverse_key)
                    if rxn_specs[reverse_key]["cat"] == cat
                        continue
                    end
                end

                # only used for mutation, to ensure that the randomly newly generated reaction is not a reaction that was previously in the reaction network
                if !isempty(delete_rxns)
                    if haskey(delete_rxns, forward_key)
                        if delete_rxns[forward_key] == cat
                            continue
                        end
                    end
                end

                if !haskey(addition_specs, forward_key) && !haskey(rxn_specs, forward_key)
                    conflicting_rxns = false
                    for r in rct
                        push!(reacting_species, r)
                    end
                    for p in prd
                        push!(reacting_species, p)
                    end
                    for c in cat
                        push!(reacting_species, c)
                    end
                    nGeneratedRxns += 1
                    addition_specs[forward_key] = Dict{String, Any}("rcts" => rct, "prds" => prd, "rxn_mech" => rxn_mech, "cat" => cat)
                end
            end
        end
        if length(reacting_species) == length(species)
            missing_species = false
        else
            addition_specs = Dict()
            reacting_species = copy(reacting_species_init)
        end
    end
    return merge(addition_specs, rxn_specs)
end

function generate_model(rxn_specs::Dict{})
    variables = Vector{String}()
    add_rxns = []

    rxn_counter = 0
    for (key, spec) in rxn_specs
        rxn_num = "J$rxn_counter"
        rxn_mech = spec["rxn_mech"]
        rcts = spec["rcts"]
        prds = spec["prds"]
        cat = spec["cat"]

        rl, vars = generate_rate_law(rcts, prds, cat, rxn_mech, rxn_counter)
        append!(variables, vars)
        add_rxn = [rxn_num, rcts, prds, rl]
        println(typeof(add_rxn))
        push!(add_rxns, add_rxn)
        rxn_counter += 1
    end
    println("final: ", variables)
    return [add_rxns, variables]
end

function generate_rate_law(rcts, prds, cat, rxn_mech::String, rxn_counter::Int64)
    add_rxn = String[]
    if rxn_mech == "UNIUNI"
        rate_law, variables = michealis_menten_rate_law(cat[1], rcts[1], prds[1], rxn_counter)
    else
        variables = String[]
        forward_rate_const = "k0_J$rxn_counter"
        reverse_rate_const = "k1_J$rxn_counter"
        push!(variables, forward_rate_const)
        push!(variables, reverse_rate_const)
        rate_law = "$forward_rate_const * "
        for (idx, rct) in enumerate(rcts)
            rate_law = string(rate_law, rct)
            if idx < length(rcts)
                rate_law = string(rate_law, " * ")
            end
        end
        rate_law = string(rate_law, " - $reverse_rate_const * ")
        for (idx, prd) in enumerate(prds)
            rate_law = string(rate_law, prd)
            if idx < length(prds)
                rate_law = string(rate_law, " * ")
            end
        end
    end
    return [rate_law, variables]
end

function michealis_menten_rate_law(E::String, S::String, P::String, rxn_counter::Int64)
    K_m = string("Km_J", rxn_counter)
    k_cat = string("kcat_J", rxn_counter)
    rate_law = "($k_cat * $E * $S) / ($S + $K_m)"
    variables = [K_m, k_cat]
    return [rate_law, variables]
end

function convert_rctprd_to_keys(rcts::Vector{String}, prds::Vector{String})
    function buildKey(species) # helper function
        species_copy = sort(species)
        key = ""
        for (idx, s) in enumerate(species_copy)
            key *= string(s)
            if idx < length(species_copy)
                key *= "_"
            end
        end
        return key
    end

    if length(rcts) == 0 || length(prds) == 0
        error("empty reactant or product list")
    end
    rct_key = buildKey(rcts)
    prd_key = buildKey(prds)
    forward_key = string(rct_key, ", ", prd_key)
    reverse_key = string(prd_key, ", ", rct_key)
    return [forward_key, reverse_key]
end

function find_boundary_species(rxn_specs::Dict{})
    edge_counter = Dict{String, Vector{Int64}}()
    boundary_species = String[]
    for (key, spec) in rxn_specs
        rcts = spec["rcts"]
        prds = spec["prds"]
        cat = spec["cat"]
        if length(cat) > 0
            for c in cat
                if !haskey(edge_counter, c)
                    edge_counter[c] = [0, 0]
                end
            end
        end

        for rct in rcts
            if haskey(edge_counter, rct)
                edge_counter[rct][1] += 1
            else
                edge_counter[rct] = [1, 0]
            end
        end
        for prd in prds
            if haskey(edge_counter, prd)
                edge_counter[prd][2] += 1
            else
                edge_counter[prd] = [0, 1]
            end
        end
    end
    for (edge, counter) in edge_counter
        if counter[1] == 0 || counter[2] == 0
            push!(boundary_species, edge)
        end
    end
    return boundary_species
end

function replace_rxn(rxn_specs::Dict{}, mutate_num::Int64, species::Vector{String}, gene_species::Vector{String}, nRxns::Int64)
    rxn_specs_copy = deepcopy(rxn_specs)
    reacting_species = Set{String}()
    delete_keys = sample(collect(keys(rxn_specs)), mutate_num, replace=false)
    delete_rxns = Dict()
    for k in delete_keys
        delete_rxns[k] = rxn_specs_copy[k]["cat"]
        delete!(rxn_specs_copy, k)
    end
    for (key, rxn) in rxn_specs_copy
        rcts = rxn["rcts"]
        prds = rxn["prds"]
        cat = rxn["cat"]
        for rct in rcts
            push!(reacting_species, rct)
        end
        for prd in prds
            push!(reacting_species, prd)
        end
        for c in cat
            push!(reacting_species, c)
        end
    end
    mutated_specs = randomize_rxns(species, gene_species, mutate_num, rxn_specs = rxn_specs_copy, reacting_species = reacting_species, delete_rxns = delete_rxns)
    return [mutated_specs, delete_rxns]
end
end
# function setBoundary(rr::Ptr{Nothing}, sid::String, boundaryCondition::Bool, forceRegen::Bool)
#   status = false
#   if forceRegen == true
#     status = ccall(dlsym(rrlib, :setBoundary), cdecl, Bool, (Ptr{Nothing}, Ptr{UInt8}, Bool), rr, sid, boundaryCondition)
#   else
#     status = ccall(dlsym(rrlib, :setBoundary), cdecl, Bool, (Ptr{Nothing}, Ptr{UInt8}, Bool), rr, pid, boundaryCondition)
#   end
#   if status == false
#     error(getLastError())
#   end
# end
