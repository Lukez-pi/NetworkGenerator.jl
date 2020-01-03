include("./NetworkGenerator.jl")
import .NetworkGenerator
include("./JuliaRoadRunnerAPI.jl")

species = ["S$i" for i = 1:7]
gene_species = species[1:2]
rxn_specs = NetworkGenerator.randomize_rxns(species, gene_species, 8)
NetworkGenerator.print_dict(rxn_specs)


add_rxns, variables = NetworkGenerator.generate_model(rxn_specs, true)

sbmlFile = "C:/Users/lukez/OneDrive/Desktop/Network Generator/Network-Generator/minimal.xml"
f = open(sbmlFile)
sbmlStr = read(f,String)
close(f)

rr = createRRInstance()         # Start up roadRunner
loadSBML(rr, sbmlStr)
addCompartment(rr, "compartment", 1.0, false)

for var in variables
    addParameter(rr, var, 0.1, false)
    println(var)
end

for s in species
    addSpecies(rr, s, "compartment", 0.1, "concentration", false)
end

for (idx, rxn) in enumerate(add_rxns)
    rid = rxn[1]
    reactants = rxn[2]
    products = rxn[3]
    kinetic_law = rxn[4]
    if length(reactants) == 1 && length(products) == 1
        temp = split(kinetic_law, "; ")
        intermediate = String(temp[1])
        rxn1 = String(temp[2])
        println("this is rxn1: ", rxn1)
        rxn2 = String(temp[3])
        println("this is intermediate: $(intermediate)")
        addSpecies(rr, intermediate, "compartment", 0.0, "concentration", false)
        addReaction(rr, string(rid, "_1"), reactants, [intermediate], rxn1, false)
        addReaction(rr, string(rid, "_2"), [intermediate], products, rxn2, false)
    else
        addReaction(rr, rid, reactants, products, kinetic_law, false)
    end
end

println("before adding boundary species: ", getNumberOfBoundarySpecies(rr))
boundary_species = NetworkGenerator.find_boundary_species(rxn_specs)


for (idx, b) in enumerate(boundary_species)
    println("boundary species: ", b)
    idx == length(boundary_species) ? regen = true : regen = false
    setBoundary(rr, b, true, regen)
end
println("after adding boundary species: ", getNumberOfBoundarySpecies(rr))
println("these are the floating species: ", getFloatingSpeciesIds(rr))
