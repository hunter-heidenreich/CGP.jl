# The various exports to be included in the 
export Chromosome,          # The abstract Chromosome type
    process,                # Processes a Chromosome
    find_active,            # Gets the list of active Nodes
    get_genes,              # Gets the genes from a Chromosome
    clone,                  # Clones a Chromrosome
    forward_connections,    # Gets a list of the forward connections
    get_output_trace        # Gets an output trace of the Chromrosome

# An abstract chromosome, to be defined
abstract type Chromosome end

#  An abstract node to be defined
abstract type Node end

# What is a CGP node look like?
mutable struct CGPNode <: Node
    connections::Array{Int64} # An array of connections coming in
    f::Function # The explicit function for the Node
    active::Bool # Whether or not it is currently active
    p::Float64 # A value to parameterize the node and connectivity
    output::Any # What output value it gives, if any
end

# Series of constuctors for CGP Nodes

# Basic constrcutor, giving inputs and a function
function CGPNode(ins::Array{Int64}, f::Function)
    CGPNode(ins, f, false, 1.0, 0.0)
end

# Extension, allowing it to be defined as active or not
function CGPNode(ins::Array{Int64}, f::Function, active::Bool)
    CGPNode(ins, f, active, 1.0, 0.0)
end

# Extension, adding a parameter value to be added
function CGPNode(ins::Array{Int64}, f::Function, active::Bool, p::Float64)
    CGPNode(ins, f, active, p, 0.0)
end

# I think this function snaps all values into the proper input range
# and then processes them
function snap(fc::Array{Float64}, p::Array{Float64})::Array{Int64}
    map(x->indmin(abs.(p-x)), fc)
end

# Processing a chromosome (individual ?) and inputs
function process(c::Chromosome, inps::Array)::Array{Float64}
    # For all input Nodes
    for i in 1:c.nin
        # Because we know these are all input values, we set their output to
        # the input value in the full CGP graph
        c.nodes[i].output = inps[i]
    end

    # Set the max length to the maximum length any given input value
    maxlength = maximum(map(length, inps))

    # Now we need to evaluate all the Nodes
    for n in c.nodes

        # Check if the Node is an active one (else we can ignore it)
        if n.active
            # Grab the output from that active Node
            output = n.output

            # If a Node is an input Node (?)
            if n.f == Config.f_input
                # The output is some sort of indexing into it (?)
                output = Config.index_in(inps, (n.p+1.0)/2.0)
            else
                # Grab the scaling factor, set it if it doesn't exist
                p = n.p
                if ~CGP.Config.weights
                    p = 1.0
                end

                # Compute the output based on function and connections and parameter
                output = CGP.Config.scaled(p * n.f(c.nodes[n.connections[1]].output,
                                                     c.nodes[n.connections[2]].output,
                                                     n.p))
            end

            # If output doesn't exist, set it to 0
            if length(output) == 0
                output = 0.0
            # Else set it to the max length of what we can handle
            elseif length(output) > maxlength
                output = output[1:maxlength]
            end

            # Update the Node's output
            n.output = output
        end
    end

    # Convert the outputs to an array of outputs
    outs = Array{Float64}(c.nout)

    # Because outputs have to be singular values, reduce to the mean
    for i in eachindex(outs)
        outs[i] = mean(c.nodes[c.outputs[i]].output)
    end

    # Return them out
    outs
end

# Function for activating Nodes with recurrences
function recur_active!(active::BitArray, connections::Array{Int64}, ind::Int64)::Void
    # Only if that Node is not active
    if ~active[ind]
        # Activate it
        active[ind] = true

        # Look at its connections
        for i in 1:2
            recur_active!(active, connections, connections[i, ind])
        end
    end
end

# Find the active sub-graph
function find_active(nin::Int64, outputs::Array{Int64}, connections::Array{Int64})::BitArray
    # Set everything to false
    active = falses(size(connections, 2)+nin)

    # Activate all the input Nodes
    active[1:nin] = true

    # For each of the final outputs, activate its connections
    for i in eachindex(outputs)
        recur_active!(active, connections, outputs[i])
    end

    # Set the inputs to false because we don't need that anymore
    active[1:nin] = false

    # Return them out
    active
end

# Another function for getting active subgraph
# Same as above, just doesn't deal with inputs
function find_active(outputs::Array{Int64}, connections::Array{Int64})::BitArray
    active = falses(size(connections, 2))
    for i in eachindex(outputs)
        recur_active!(active, connections, outputs[i])
    end
    active
end

# Resets a chromrosome
function reset!(c::Chromosome)
    for n in c.nodes
        n.output = 0.0
    end
end

# Has to implementeed down the line
function get_positions(c::Chromosome)
    error("Must be implemented in subclass")
end

# A recursive function for tracing output
function recur_output_trace(c::Chromosome, ind::Int64, visited::Array{Int64})
    # Check if a Node has been visited
    if ~contains(==, visited, ind)

        # add it to the list of visited Nodes
        append!(visited, [ind])

        # If the Node is not an input Node and is not an Input function
        if (ind > c.nin) && (c.nodes[ind].f != Config.f_input)

            # Visit all its connections
            for i in c.nodes[ind].connections
                recur_output_trace(c, i, visited)
            end
        end
    end

    # Return out the visited Nodes
    visited
end

# Calls the output trace for a given chromosome and a specific output index
function get_output_trace(c::Chromosome, output_ind::Int64)
    # similar to decode, just return a list of node indices that determine the output
    recur_output_trace(c, c.outputs[output_ind], Array{Int64}(0))
end

# Gets a full output trace
function get_output_trace(c::Chromosome, outputs::Array{Int64})
    if length(outputs) > 0
        return unique(reduce(vcat, map(x->get_output_trace(c, x), outputs)))
    else
        return Array{Int64}(0)
    end
end

# Gets a full output trace by implicitly collecting its outputs
function get_output_trace(c::Chromosome)
    # same as indexing over active nodes
    get_output_trace(c, collect(1:c.nout))
end

# Is a subclass implementation
function node_genes(c::Chromosome)
    # number of genes per node (position, x, y, f, param)
    error("Must be implemented in subclass")
end

# Gets all the genes for a chromosome and Node
function get_genes(c::Chromosome, node_id::Int64)
    c.genes[(c.nin+c.nout)+((node_id-1-c.nin)*node_genes(c))+(1:node_genes(c))]
end

# Gets all the genes for a chromosome and an array of Nodes
function get_genes(c::Chromosome, nodes::Array{Int64})
    if length(nodes) > 0
        return reduce(vcat, map(x->get_genes(c, x), nodes))
    else
        return Array{Int64}(0)
    end
end

# Updates Genes for a chromosome at a given Node
function set_genes!(c::Chromosome, node_id::Int64, genes::Array{Float64})
    c.genes[(c.nin+c.nout)+((node_id-1-c.nin)*node_genes(c))+(1:node_genes(c))] = genes
end

# Gets all the forward connections
function forward_connections(c::Chromosome)
    # Get the list of connections for all Nodes
    connections = [[i] for i in 1:length(c.nodes)]

    #  Now iterate over each Node
    for ci in eachindex(c.nodes)
        # Iterate over its connections
        for i in 1:2
            # Get the connection
            conn = c.nodes[ci].connections[i]

            # IF it's greater than 0, we want it (it exists)
            if conn > 0
                # Check if it's already in the connection list
                # If no, include it
                if ~(contains(==, connections[ci], conn))
                    append!(connections[ci], [conn])
                end

                # For each additional connection, check if its in the list
                # and include it
                for j in eachindex(connections[conn])
                    if ~(contains(==, connections[ci], j))
                        append!(connections[ci], [j])
                    end
                end
            end
        end
    end

    # Return out all the connections
    connections
end

# Clones a chromosome
function clone(c::Chromosome)
    # call gene constructor method
    typeof(c)(deepcopy(c.genes), c.nin, c.nout)
end
