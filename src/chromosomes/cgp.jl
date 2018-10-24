export CGPChromo, node_genes, get_positions

# vanilla CGP

# Inherits from base Chromosome
struct CGPChromo <: Chromosome
    genes::Array{Float64}
    nodes::Array{CGPNode}
    outputs::Array{Int64}
    nin::Int64
    nout::Int64
end

# Constructor 1: Give the genes, inputs, and outputs
function CGPChromo(genes::Array{Float64}, nin::Int64, nout::Int64)::CGPChromo
    # Comput the number of nodes
    num_nodes = Int64(ceil((length(genes)-nin-nout)/4))

    # These are genes that don't evolve (?)
    rgenes = reshape(genes[(nin+nout+1):end], (4, num_nodes))'

    #  Compute all the connections
    connections = Array{Int64}(2, nin+num_nodes)
    connections[:, 1:nin] = zeros(2, nin)

    #  Compute the positions
    positions = collect(1:(nin+num_nodes))/(1.0*(nin+num_nodes))

    # Not sure what's happening in this block
    fc = deepcopy(hcat(zeros(2, nin), [rgenes[:, 2]'; rgenes[:, 3]']))
    e = (1.0 .- positions) .* Config.recurrency .+ positions
    for j in 1:size(fc)[1]
        fc[j, :] .*= e
    end

    # If no recurrency, make sure of that
    if Config.recurrency == 0
        for i in (nin+1):length(positions)
            for j in 1:size(fc)[1]
                fc[j, i] = max(fc[j,i], positions[i-1])
            end
        end
    end

    # Snap all the connections in place
    connections = snap(fc, positions)

    # Set all the functions
    functions = Array{Function}(nin+num_nodes)
    functions[1:nin] = Config.f_input
    functions[(nin+1):end] = map(i->Config.index_in(Config.functions, i), rgenes[:, 3])

    # Set up the Outputs
    outputs = Int64.(ceil.(genes[nin+(1:nout)]*(nin+num_nodes)))

    # Grab the active graph
    active = find_active(nin, outputs, connections)

    # Compute the parameters
    params = [zeros(nin); 2.0*rgenes[:, 4]-1.0]

    # Get all teh nodes
    nodes = Array{CGPNode}(nin+num_nodes)

    # Initialize all the Nodes
    for i in 1:(nin+num_nodes)
        nodes[i] = CGPNode(connections[:, i], functions[i], active[i], params[i])
    end

    # Return out the Chromosome
    CGPChromo(genes, nodes, outputs, nin, nout)
end

# Constructor 2: Just give the inputs and outputs
function CGPChromo(nin::Int64, nout::Int64)::CGPChromo
    n_nodes = Config.static_node_size
    if Config.bloat()
        n_nodes = Config.starting_nodes
    end
    CGPChromo(rand(nin+nout+4*n_nodes), nin, nout)
end

# Constructor 3: Give a full Chromosome
function CGPChromo(c::CGPChromo)::CGPChromo
    gene_mutate(c)
end

# Returns number of genes
function node_genes(c::CGPChromo)
    4
end

# Gets the positions of the Nodes
function get_positions(c::CGPChromo)
    collect(1:length(c.nodes))/(1.0*length(c.nodes))
end

# Adds a subtree to Chromosome
function add_subtree(c::CGPChromo)
    add_nodes(c)
end

#  Deletes a subtree from the Chromosome
function delete_subtree(c::CGPChromo)
    delete_nodes(c)
end
