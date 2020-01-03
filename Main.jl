# include("./NetworkGenerator.jl")
# using NetworkGenerator
include("./JuliaRoadRunnerAPI.jl")

species = ["S$i" for i = 1:7]
gene_species = Set(["1", "2"])
rxn_specs = randomize_rxns(species, gene_species, 8)
print_dict(rxn_specs)


add_rxns, variables = generate_model(rxn_specs)
println(add_rxns)
println(variables)

sbmlFile = "C:/Users/lukez/OneDrive/Desktop/Network Generator/Network-Generator/minimal.xml"
f = open(sbmlFile)
sbmlStr = read(f,String)
close(f)

rr = createRRInstance()         # Start up roadRunner
loadSBML(rr, sbmlStr)
addCompartment(rr, "compartment", 1.0, false)

for var in variables
    addParameter(rr, var, 0.1, false)
end

for s in species
    println(typeof(s))
    println(s)
    addSpecies(rr, s, "compartment", 0.1, "concentration", false)
end

for (idx, rxn) in enumerate(add_rxns)
    rid = rxn[1]
    reactants = rxn[2]
    products = rxn[3]
    kineticLaw = rxn[4]
    addReaction(rr, rid, reactants, products, kineticLaw, false)
end

println("before adding boundary species: ", getNumberOfBoundarySpecies(rr))
boundary_species = find_boundary_species(rxn_specs)
for (idx, b) in enumerate(boundary_species)
    println(b)
    idx == length(boundary_species) ? regen = true : regen = false
    setBoundary(rr, b, true, regen)
end
println("after adding boundary species: ", getNumberOfBoundarySpecies(rr))
