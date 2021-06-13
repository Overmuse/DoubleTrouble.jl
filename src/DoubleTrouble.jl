module DoubleTrouble

using CSV
using DataFrames
using Dates
using Impute: Impute, impute, Interpolate
using JSON3
using Optim: BFGS, Fminbox, optimize, minimizer
using Roots: find_zero
using ShiftedArrays: ShiftedArray
using Statistics: quantile, var
using StatsBase: ordinalrank
using TimeZones: ZonedDateTime, astimezone, @tz_str

export main

include("utils.jl")
include("data.jl")
include("parameter_tuning.jl")
include("selection.jl")

function main()
    @info "Starting DoubleTrouble.jl"
    data = read_adjusted_prices(ENV["IN_FILE"])
    tickers = names(data)
    shortterm_spreads = spreads(data[end-30*79:end, :])
    longterm_data = data[reduce(vcat, [[79*(j-1)+1, 79j] for j in 1:100]), :]
    longterm_spreads = spreads(longterm_data)
    initial_prices = data[1, :]
    pairs_data = generate_trade_pairs(tickers, shortterm_spreads, longterm_spreads, initial_prices)
    CSV.write(ENV["OUT_FILE"], DataFrame(pairs_data))
end

end # module
