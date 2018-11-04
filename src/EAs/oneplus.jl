# The file for defining what a one plus evolution scheme looks like

export oneplus

# Pass it your inputs, outputs, and fitness function
# The default Chromosome type is CGP, the seed is fixed, and there should be no
# expert
# String ID is optional
function oneplus(nin::Int64, nout::Int64, fitness::Function;
                 ctype::DataType=CGPChromo, seed::Int64=0, expert::Any=nothing,
                 id::String="")

    # Set up a population of CGP individuals
    population = Array{ctype}(Config.lambda)

    # For each of them, set them up with a Chromosome
    for i in eachindex(population)
        # This line actually creates random individuals
        population[i] = ctype(nin, nout)
    end

    # If we are using an expert, plug it in here
    if expert != nothing
        population[1] = expert
    end

    # We assume the best one is the first one
    best = population[1]

    # Setting some basic parameters
    max_fit = -Inf
    eval_count = 0
    fits = -Inf*ones(Float64, Config.lambda)

    # While we have less than the number of total evaluations we are going to
    # do
    while eval_count < Config.total_evals
        # evaluation
        log_gen = false

        # Loop through each individual
        for p in eachindex(population)
            # If it doesn't have a fitness set yet
            if fits[p] == -Inf
                # Evaluate its fitness
                fit = fitness(population[p])
                eval_count += 1

                # If it is the best or equal to the best, we clone it and save it
                if fit >= max_fit
                    best = clone(population[p])
                    if fit > max_fit
                        max_fit = fit

                        # We only log generations where candidates exceed
                        # previous performance
                        log_gen = true
                    end
                end

                # Update the fitness
                fits[p] = fit

                # Log if that was the last generation
                if eval_count == Config.total_evals
                    log_gen = true
                    break
                end
            end
        end

        # Evaluate that generation
        # eval(Config.log_function)(id, seed, eval_count, max_fit, best, fitness, ctype, log_gen)

        if eval_count == Config.total_evals
            break
        end

        # selection
        fits .= -Inf
        for p in eachindex(population)
            population[p] = mutate(best)
        end

        # size limit
        for i in eachindex(population)
            # If a mutated individual exceed what we expect,
            # Generate a random one (?)
            # Side Note: Shouldn't we just generate one off the best?
            if length(population[i].nodes) > Config.node_size_cap
                population[i] = ctype(nin, nout)
                fits[i] = -Inf
            end
        end
    end

    # When done, return out the max fit candidate and its Chromosome
    max_fit, best.genes
end
