function read_adjusted_prices(file)
    str = read(file, String)
    dateformat = dateformat"yyyy-mm-ddTHH:MM:SS\Z"
    parsed = JSON3.read(str, Dict{String, Tuple{Vector{DateTime}, Vector{Float64}}}; dateformat)
    data = DataFrame()
    for (ticker, (dates, prices)) in parsed
        temp = DataFrame(ticker = ticker, datetime = astimezone.(ZonedDateTime.(dates, tz"UTC"), tz"America/New_York"), price = prices)
        data = vcat(data, temp)
    end
    data = sort(unstack(data, :ticker, :price), :datetime)
    filter!(data) do row
        (Time(row.datetime) >= Time(9, 30)) && (Time(row.datetime) <= Time(16, 00))
    end
    data.date = Date.(data.datetime)
    data = combine(groupby(data, :date)) do group
        mapcols(group) do col
            if all(ismissing, col)
                col
            else
                locf(impute(col, Interpolate(; limit = 1)))
            end
        end
    end
    disallowmissing(Impute.filter(data, dims = :cols)[:, Not([:datetime, :date])])
end

function get_open_prices(tickers)
    creds = get_credentials()
    map(tickers) do ticker
        Float64(get_snapshot(creds, ticker)["day"]["o"])
    end
end
