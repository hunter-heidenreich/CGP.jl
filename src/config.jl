# This is the Config Module
module Config

# All them imports
using YAML
using Logging
using PaddedViews
using Distributions
using ArgParse

include("functions.jl")
functions = Array{Function}(0)
# All them imports

# Initialize a configuration based on a read dictionary
function init(config::Dict)
    # Look through all the keys in the dictionary
    for k in keys(config)
        # If it's a function list, add it and load the function
        if k == "functions"
            append!(functions, load_functions(config["functions"]))
        else
            # If it's not Nothing, then parse it directly
            if config[k] != nothing
                eval(parse(string(k, "=", config[k])))
            end
        end
    end
end

# Looks like this is some sort of or operation for checking things?
function bloat()
    ((mutate_method in [:mixed_node_mutate, :adaptive_node_mutate,
                        :mixed_subtree_mutate, :adaptive_subtree_mutate]) ||
     (crossover_method in [:output_graph_crossover, :subgraph_crossover]))
end

# Load a config file
function init(file::String)
    init(YAML.load_file(file))
end

# Resest all the available functions
function reset()
    empty!(functions)
end

# Add the settings
function add_arg_settings!(s::ArgParseSettings)

    mutations = [":gene_mutate", ":mixed_node_mutate", ":mixed_subtree_mutate",
                 ":adaptive_node_mutate", ":adaptive_subtree_mutate"]

    crossovers = [":single_point_crossover", ":random_node_crossover",
                  ":aligned_node_crossover", ":proportional_crossover",
                  ":output_graph_crossover", ":subgraph_crossover"]

    distances = [":genetic_distance", ":positional_distance",
                 ":constant_functional_distance", ":random_functional_distance",
                 ":active_distance"]

    # Restrict various parameters to what we have implemented
    @add_arg_table s begin
        "--mutate_method"
            default = nothing
            range_tester = (x->x ∈ mutations)
            help = "mutation method; must be one of " * join(mutations, ", ", " or ")
        "--crossover_method"
            default = nothing
            range_tester = (x->x ∈ crossovers)
            help = "crossover method; must be one of " * join(crossovers, ", ", " or ")
        "--distance_method"
            default = nothing
            range_tester = (x->x ∈ distances)
            help = "distance method; must be one of " * join(distances, ", ", " or ")
    end

    # Other parameters that will be encoded
    params = ["input_start", "recurrency", "input_mutation_rate",
        "output_mutation_rate", "node_mutation_rate", "node_size_delta",
        "modify_mutation_rate", "ga_elitism_rate", "ga_crossover_rate",
        "ga_mutation_rate"]

    # Load up the float parameters
    for p in params
        add_arg_table(s, ["--$p"], Dict(:help=>"Parameter: $p", :arg_type=>Float64))
    end

    # Load up the int parameters
    for p in ["lambda", "ga_population", "starting_nodes", "static_node_size",
              "node_size_cap", "total_evals"]
        add_arg_table(s, ["--$p"], Dict(:help=>"Parameter: $p", :arg_type=>Int64))
    end

    # Load up the bool parameters
    for p in ["active_mutate", "weights", "save_best"]
        add_arg_table(s, ["--$p"], Dict(:help=>"Parameter: $p", :arg_type=>Bool))
    end

    # Return them
    s
end

# Get the settings and set them up
function get_arg_settings()
    s = ArgParseSettings()
    add_arg_settings!(s)
    s
end

# Print all the settings to a string
function to_string()
    @sprintf(
        "%s %s %s %s %d %d %d %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %d %d %0.3f %0.3f %0.3f",
        string(mutate_method), string(active_mutate), string(crossover_method),
        string(weights), starting_nodes, static_node_size, node_size_cap,
        input_start, recurrency, input_mutation_rate,
        output_mutation_rate, node_mutation_rate, node_size_delta,
        modify_mutation_rate, lambda, ga_population, ga_elitism_rate,
        ga_crossover_rate, ga_mutation_rate)

end

# append!(functions, [f_input])
# All we need to actually export out of this Module is
# A function to initialize
# A function to reset
export init, reset
end
