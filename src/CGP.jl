module CGP

# This defines what the import module of 'CGP is'

using Logging
using PaddedViews
using Distributions

Logging.configure(level=INFO)

# Various includes we wish to have in this module

include("config.jl")
include("chromosome.jl")
include("distance.jl")
include("mutation.jl")
include("crossover.jl")
include("chromosomes/cgp.jl")
include("chromosomes/pcgp.jl")
include("logging.jl")
include("EAs/oneplus.jl")
include("EAs/cgpneat.jl")
include("EAs/ga.jl")
include("EAs/cmaes.jl")

# The series of evolutionary strategies to be included and used
EAs = [oneplus, GA]

#  The types of chromosomes, normal CGP and positional CGP
chromosomes = [CGPChromo, PCGPChromo]

# The types of mutation strategies that are included in the package
mutations = [:gene_mutate, :mixed_node_mutate, :mixed_subtree_mutate, :adaptive_node_mutate,
             :adaptive_subtree_mutate]

# The types of crossover strategies included in the package
crossovers = [:single_point_crossover, :random_node_crossover, :aligned_node_crossover,
              :proportional_crossover, :output_graph_crossover, :subgraph_crossover]

# The distance measures included in the package
distances = [:genetic_distance, :positional_distance, :constant_functional_distance,
             :random_functional_distance, :active_distance]

end
