using ArcadeLearningEnvironment
using CGP
using Logging
using ArgParse
import Images

include("../graphing/graph_utils.jl")

CGP.Config.init("cfg/atari.yaml")

function stop_playing_check(cur_reward::Float64, frame_step::Int64)
    global max_reward_dict

    if haskey(max_reward_dict, frame_step)
        print("[", frame_step, "] Comparing: ", cur_reward, " to ", max_reward_dict[frame_step])

        if max_reward_dict[frame_step] == 0
            println()
            return false
        end

        val = ((max_reward_dict[frame_step] - cur_reward) / max_reward_dict[frame_step])
        println(" -- ", val)

        if val > 0.5
            return true
        end
    end

    false
end

@everywhere begin
    function score_velocity_score(reward_dict::Dict)
        score_vel = 0.0
        prev = 0.0
        step_scale = 1000

        for step in reward_dict
            if step[1] > 0 && step[2] != -Inf
                cur = step[2]

                diff = cur - prev
                score_vel += diff / step_scale

                prev = cur
            end
        end

        score_vel
    end

    function play_atari(c::Chromosome, id::String, seed::Int64;
                        render::Bool=false, folder::String=".", max_frames=18000,
                        step_size=200)

        score_dict = Dict()
        # Get the game with ID and seed
        game = Game(id, seed)

        # Reset the seed
        seed_reset = rand(0:100000)
        srand(seed)

        # Reset all parameters
        reward = 0.0
        frames = 0
        p_action = game.actions[1]
        outputs = zeros(Int64, c.nout)

        # Keep playing, while the game isn't over
        while ~game_over(game.ale)

            # Have the Chromosome process the current inputs
            output = process(c, get_rgb(game))

            # Log the move count
            outputs[indmax(output)] += 1

            # Do that action
            action = game.actions[indmax(output)]
            reward += act(game.ale, action)

            # Possibly repeat action randomly
            if rand() < 0.25
                reward += act(game.ale, action) # repeat action here for seeding
            end

            # If we are rendering the game, save the screen to a file
            if render
                screen = draw(game)
                filename = string(folder, "/", @sprintf("frame_%06d.png", frames))
                Images.save(filename, screen)
            end

            if (frames % step_size) == 0
                score_dict[frames] = reward
            end

            # Always log the frames
            frames += 1
            if frames > max_frames
                Logging.debug(string("Termination due to frame count on ", id))
                break
            end
        end

        # Close the game and reset the seed
        close!(game)
        srand(seed_reset)

        score_dict[frames] = reward

        vel_score = score_velocity_score(score_dict)

        # Return the reward and the list of outputs
        reward, outputs, vel_score
    end
end

# Parses all the command line arguments
function get_args()
    s = ArgParseSettings()

    @add_arg_table(
        s,
        "--seed", arg_type = Int, default = 0,
        "--log", arg_type = String, default = "atari.log",
        "--id", arg_type = String, default = "centipede",
        "--ea", arg_type = String, default = "oneplus",
        "--chromosome", arg_type = String, default = "CGPChromo",
        "--frames", arg_type = Int, default = 18000,
        "--render", action = :store_const, constant = true, default = false,
    )

    parse_args(CGP.Config.add_arg_settings!(s))
end

# Gets the best Chromrosomes from a log file
function get_bests(logfile::String)
    bests = []
    for line in readlines(logfile)
        if contains(line, ":C:")
            genes = eval(parse(split(line, ":C: ")[2]))
            append!(bests, [genes])
        end
    end
    bests
end

# Gets the input/output parameters
function get_params(args::Dict)
    game = Game(args["id"], args["seed"])
    nin = 3 # r g b
    nout = length(game.actions)
    close!(game)
    nin, nout
end

# Renders a gene
function render_genes(genes::Array{Float64}, args::Dict;
                    ctype::DataType=CGPChromo, id::Int64=0)

    # Gets the Chromrosome
    nin, nout = get_params(args)
    chromo = ctype(genes, nin, nout)

    # Selects the folder
    folder = string("frames/", args["id"], "_", args["seed"], "_", id)
    mkpath(folder)

    # Plays the Atari with that individual
    reward, out_counts, fit_score = play_atari(chromo, args["id"], args["seed"];
                                        render=false, folder=folder,
                                        max_frames=args["frames"])

    # Logs that individual's performance
    println(@sprintf("R: %s %d %d %0.5f %d %d",
                     args["id"], args["seed"], id, reward,
                     sum([n.active for n in chromo.nodes]),
                     length(chromo.nodes)))

    Logging.info(string("Top individual scores: "))
    Logging.info(string("\tFitness Score: ", fit_score))
    Logging.info(string("\tGame Score: ", reward))

    println("Top individual scores: ")
    println("\tFitness Score: ", fit_score)
    println("\tGame Score: ", reward)

    # Copy our genes
    new_genes = deepcopy(chromo.genes)

    # Not sure what this all is about
    for o in 1:nout
        if out_counts[o] == 0.0
            # disable this section of the graph by forcing
            # the output to connect to an input node
            new_genes[nin+o] = 0.00001
        end
    end

    # Check how many active outputs we actually had
    active_outputs = out_counts .> 0

    # Take the new Chromosome with slight modification to draw
    chromo2 = ctype(new_genes, nin, nout)
    file =  string("graphs/", args["id"], "_", args["seed"], "_", id, ".pdf");
    chromo_draw(chromo2, file, active_outputs=active_outputs)
    # chromo_draw(chromo2)
end

# We run this as a non-interactive experiment
if ~isinteractive()

    # Get the args
    args = get_args()

    # Configure the CGP experiment
    CGP.Config.init(Dict([k=>args[k] for k in setdiff(
        keys(args), ["seed", "log", "id", "ea", "chromosome"])]...))

    # Seed the random generator
    srand(args["seed"])

    # Set up logging
    Logging.configure(filename=string(args["id"], ".log"), level=INFO)

    # Get the inputs/outputs
    nin, nout = get_params(args)

    # Select evolution schema
    ea = eval(parse(args["ea"]))

    # Select chromosome type
    ctype = eval(parse(args["chromosome"]))

    # Define the fitness function
    # We use [1] here because it selects the reward as the fitness
    # measure
    fit = x->play_atari(x, args["id"], args["seed"];
                            max_frames=args["frames"])[3]

    # Actually run the EA experiment
    maxfit, best = ea(nin, nout, fit;
                      seed=args["seed"], id=args["id"], ctype=ctype)

    # Log the maxfit individual
    Logging.info(@sprintf("E%0.6f", -maxfit))

    # If we are rendering, render this best Chromosome
    if args["render"]
        render_genes(best, args; ctype=ctype)
    end
end
