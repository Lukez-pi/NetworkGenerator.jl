using Test
include("./NetworkGenerator.jl")
using .NetworkGenerator

@testset "rate_law" begin
    # UNIUNI
    rl_UU, var_UU = NetworkGenerator.generate_rate_law(["S1"], ["S3"], ["E5"], "UNIUNI", 2, false)
    @test rl_UU == "(kcat_J2 * E5 * S1) / (S1 + Km_J2)"
    @test Set(var_UU) == Set(["kcat_J2", "Km_J2"])

    rl_UU, var_UU = NetworkGenerator.generate_rate_law(["S1"], ["S3"], ["E5"], "UNIUNI", 2, true)
    temp = split(rl_UU, "; ")
    @test temp[1] == "inter_E5-S1"
    @test temp[2] == "k1_J2 * S1 - k-1_J2 * $(temp[1])"
    @test temp[3] == "kcat_J2 * $(temp[1])"
    @test Set(var_UU) == Set(["k1_J2", "k-1_J2", "kcat_J2"])

    #UNIBI
    rl_UB, var_UB = NetworkGenerator.generate_rate_law(["P3"], ["K51", "5000"], [], "BIUNI", 92, false)
    @test rl_UB == "k0_J92 * P3 - k1_J92 * K51 * 5000"
    @test Set(var_UB) == Set(["k0_J92", "k1_J92"])

    rl_UB, var_UB = NetworkGenerator.generate_rate_law(["P3"], ["K51", "5000"], [], "BIUNI", 92, true)
    @test rl_UB == "k0_J92 * P3 - k1_J92 * K51 * 5000"
    @test Set(var_UB) == Set(["k0_J92", "k1_J92"])

    # BIUNI
    rl_BU, var_BU = NetworkGenerator.generate_rate_law(["S1", "S2"], ["S3"], [], "BIUNI", 2, false)
    @test rl_BU == "k0_J2 * S1 * S2 - k1_J2 * S3"
    @test Set(var_BU) == Set(["k0_J2", "k1_J2"])

    rl_BU, var_BU = NetworkGenerator.generate_rate_law(["S1", "S2"], ["S3"], [], "BIUNI", 2, true)
    @test rl_BU == "k0_J2 * S1 * S2 - k1_J2 * S3"
    @test Set(var_BU) == Set(["k0_J2", "k1_J2"])
    # BIBI
    rl_BB, var_BB = NetworkGenerator.generate_rate_law(["S7", "M2"], ["D3", "6"], [], "BIBI", 7, false)
    @test rl_BB == "k0_J7 * S7 * M2 - k1_J7 * D3 * 6"
    @test Set(var_BB) == Set(["k0_J7", "k1_J7"])

    rl_BB, var_BB = NetworkGenerator.generate_rate_law(["S7", "M2"], ["D3", "6"], [], "BIBI", 7, true)
    @test rl_BB == "k0_J7 * S7 * M2 - k1_J7 * D3 * 6"
    @test Set(var_BB) == Set(["k0_J7", "k1_J7"])

    rl_MM, var_MM = NetworkGenerator.michealis_menten_rate_law("E1", "S2", "P3", 4)
    @test rl_MM == "(kcat_J4 * E1 * S2) / (S2 + Km_J4)"
    @test Set(var_MM) == Set(["kcat_J4", "Km_J4"])
end

@testset "find_boundary_species" begin
    rxn_spec1 = Dict(
        "1_4, 5_7" => Dict("rcts" => ["4", "1"], "cat" => Any[], "rxn_mech" => "BIBI", "prds" => ["5", "7"]),
        "3, 4_7" => Dict("rcts" => ["3"], "cat" => Any[], "rxn_mech" => "UNIBI", "prds" => ["4", "7"]),
        "1_4, 2" => Dict("rcts" => ["1", "4"], "cat" => Any[], "rxn_mech" => "BIUNI", "prds" => ["2"]),
        "1_6, 5" => Dict("rcts" => ["6", "1"], "cat" => Any[], "rxn_mech" => "BIUNI", "prds" => ["5"]),
        "3_7, 2_5" => Dict("rcts" => ["3", "7"], "cat" => Any[], "rxn_mech" => "BIBI", "prds" => ["5", "2"]),
        "4, 3_7" => Dict("rcts" => ["4"], "cat" => Any[], "rxn_mech" => "UNIBI", "prds" => ["7", "3"]),
        "2_6, 1_3" => Dict("rcts" => ["6", "2"], "cat" => Any[], "rxn_mech" => "BIBI", "prds" => ["1", "3"]),
        "1_3, 5_7" => Dict("rcts" => ["3", "1"], "cat" => Any[], "rxn_mech" => "BIBI", "prds" => ["5", "7"]))
    @test Set(["5", "6"]) == Set(NetworkGenerator.find_boundary_species(rxn_spec1))

    rxn_spec2 = Dict(
        "1, 2" => Dict("rcts" => ["1"], "cat" => ["4"], "rxn_mech" => "UNIUNI", "prds" => ["2"])
    )
    @test Set(["1", "2", "4"]) == Set(NetworkGenerator.find_boundary_species(rxn_spec2))
end

@testset "randomize_reactions" begin
    for i = 1:1000
        n_species = rand(5:21)
        species = ["S$i" for i in 1:n_species]
        gene_species = species[1:rand(1:n_species-1)]
        nRxns = rand((n_species - 2):(2 * n_species))
        rxn_specs = NetworkGenerator.randomize_rxns(species, gene_species, nRxns)
        @test length(rxn_specs) == nRxns

        species_present = Set()
        for (key, rxn) in rxn_specs
            cat = rxn["cat"]
            rcts = rxn["rcts"]
            prds = rxn["prds"]
            for rct in rcts
                push!(species_present, rct)
            end
            for prd in prds
                push!(species_present, prd)
            end
            for c in cat
                push!(species_present, c)
            end

            split_key = split(key, ", ")
            reverse_key = "$(split_key[2]), $(split_key[1])"
            if haskey(rxn_specs, reverse_key)
                @test rxn_specs[reverse_key]["cat"] != cat[1]
            end
        end
        @test length(species_present) == n_species
    end
end

@testset "mutate_rxns" begin
    for i = 1:1000
        n_species = rand(5:21)
        species = ["S$i" for i in 1:n_species]
        gene_species = species[1:rand(1:n_species-1)]
        nRxns = rand((n_species - 2):(2 * n_species))
        mutate_num = rand(1:nRxns)
        rxn_specs = NetworkGenerator.randomize_rxns(species, gene_species, nRxns)
        mutate_specs, delete_rxns = NetworkGenerator.replace_rxn(rxn_specs, mutate_num, species, gene_species, nRxns)
        @test length(delete_rxns) == mutate_num

        @test length(mutate_specs) == nRxns

        species_present = Set()
        new_rxn_keys = Set()
        new_counter = 0
        for (key, rxn) in mutate_specs
            cat = rxn["cat"]
            rcts = rxn["rcts"]
            prds = rxn["prds"]
            for rct in rcts
                push!(species_present, rct)
            end
            for prd in prds
                push!(species_present, prd)
            end
            for c in cat
                push!(species_present, c)
            end
            split_key = split(key, ", ")
            reverse_key = "$(split_key[2]), $(split_key[1])"
            if haskey(mutate_specs, reverse_key)
                @test mutate_specs[reverse_key]["cat"] != cat[1]
            end

            if length(mutate_specs[key]["cat"]) == 0
                if in(key, new_rxn_keys)
                    println("bad key: ", key)
                end
                 push!(new_rxn_keys, key)
                 new_counter += 1
            else
                new_key = string(key, ", ", mutate_specs[key]["cat"][1])
                if in(new_key, new_rxn_keys)
                    println("bad key: ", new_key)
                end
                push!(new_rxn_keys, new_key)
                new_counter += 1
            end
        end
        @test length(species_present) == n_species

        old_rxn_keys = Set()
        for (key, rxn) in rxn_specs
            if length(rxn_specs[key]["cat"]) == 0
                if in(key, old_rxn_keys)
                    println("bad key: ", key)
                end
                push!(old_rxn_keys, key)
            else
                new_key = string(key, ", ", rxn_specs[key]["cat"][1])
                if in(new_key, old_rxn_keys)
                    println("bad key: ", new_key)
                end
                push!(old_rxn_keys, new_key)
            end
        end
        @test length(setdiff(new_rxn_keys, old_rxn_keys)) == length(setdiff(old_rxn_keys, new_rxn_keys))
        @test length(setdiff(new_rxn_keys, old_rxn_keys)) == mutate_num
    end
end
