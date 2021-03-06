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

    log_anyways = 200
    buffer = Config.lambda

    # While we have less than the number of total evaluations we are going to
    # do
    while eval_count < Config.total_evals
        # evaluation
        log_gen = false

        # We are going to loop through this using Threads.@threads
        for p in eachindex(population)
            println("Evaluating #", p + eval_count, "(", (100 * (p + eval_count)/Config.total_evals), "%)")

            if p + eval_count > Config.total_evals
                break
            end

            if fits[p] == -Inf
                # Get fitness score for each and update the array
                fit = fitness(population[p])
                # Update the fitness
                fits[p] = fit
            end
        end
        # This is the only thing we should do in this loop

        # Afterwards, get max and argmax of population
        cur_max = maximum(fits)
        cur_max_id = indmax(fits)

        println("Finished generation.")

        # If max > max_fit, then:
        if cur_max >= max_fit
            # # Update max_fit
            max_fit = cur_max
            # # Update best
            best = clone(population[cur_max_id])
            # # Flip log_gen
            log_gen = true
            println("New Max Fit: ", max_fit)
        end

        # eval_count += len(population)
        eval_count += length(population)

        # if eval_count > Config.total_evals, then
        if eval_count > Config.total_evals
            # # Flip log_gen
            log_gen = true
        end

        if eval_count % log_anyways < buffer
            log_gen = true
        end

        # Call eval
        # Evaluate that generation
        eval(Config.log_function)(id, seed, eval_count, max_fit, best, fitness, ctype, true)

        # if eval_count > Config.total_evals, then
        if eval_count >= Config.total_evals
            # # break
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
