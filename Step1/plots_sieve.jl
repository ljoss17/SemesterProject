using Plots

function plot_sieve_totality(e_t)
    plot_array = []

    for g in 200:200
        local x = []
        local y = []
        for e in e_t+25:1:e_t+40
            push!(x, e)
            push!(y, sieve_total_validity(g, e, e_t))
        end
        display(
            plot(
                x,
                y,
                #yscale = :log10,
                title=string("G : ", g, " E_T : ", e_t),
                xlabel="E",
                ylabel="ϵ-total validity for Sieve",
                legend=false,
            )
        )
    end
end
