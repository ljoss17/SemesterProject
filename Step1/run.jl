include("external_libs.jl")
include("optimize_contagion.jl")

function get_optimal_parameters(N::Int64=1024, f::Float64=0.1, bound::Float64=1e-10)
    optimize_contagion(N, f, bound)
    f = CSV.File(pwd()*"/Parameters/Contagion/contagion_params_N($N)_bound($bound).csv")
    println("G : $(f.G[1]), E : $(f.E[1]), E_thr : $(f.E_thr[1]), R : $(f.R[1]), R_thr : $(f.R_thr[1]), D : $(f.D[1]), D_thr : $(f.D_thr[1])")
end

function avg_vs_security(N::Int64=1024, f::Float64=0.1)
    bound::Float64 = 0.0
    best_ϵ::Float128 = 0.0
    x::Array{Int64, 1} = zeros(11)
    y::Array{Float128, 1} = zeros(11)
    id = 1
    for pw in -16:-6
        println("pw : $pw")
        bound = 10.0^pw
        best_ϵ, x[id] = optimize_contagion(N, f, bound)
        y[id] = best_ϵ
        id = id+1
    end
    plot(
        x,
        y,
        yscale=:log10,
        title=string("Security vs. Average sample size"),
        xlabel="Average sample size (S)",
        ylabel="Security (ϵ)",
        labels="N=$N, f=$f"
    )
    savefig(pwd()*"/avg_vs_sec_N($N)_f($f).png")
end

if length(ARGS) == 0
    println("Arguments required to run the script. Following format : run.jl <command> <arg1> <arg2> ...\nPossible commands :")
    println("\tFor average vs security : avg_vs_sec N f")
    println("\tFor optimized parameters : opt N f bound")
    exit()
end
if ARGS[1] == "avg_vs_sec"
    avg_vs_security(parse(Int64, ARGS[2]), parse(Float64, ARGS[3]))
elseif ARGS[1] == "opt"
    get_optimal_parameters(parse(Int64, ARGS[2]), parse(Float64, ARGS[3]), parse(Float64, ARGS[4]))
else
    println("Command $(ARGS[1]) doesn't exist.\nExisting commands :")
    println("\tFor average vs security : avg_vs_sec N f")
    println("\tFor optimized parameters : opt N f bound")
    exit()
end
