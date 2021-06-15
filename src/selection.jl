function score(θ, σ, δ₁)
    σ / (2θ) * (1 - exp(-2θ * δ₁/78))
end

function L_score(θ, σ)
    θ < 0 && return Inf
    σ / (2θ) * (1 - exp(-2θ*(1/252)))
end

function select_pairs(L_scores, S_scores, n=50)
    ranking = ordinalrank(L_scores) .+ ordinalrank(S_scores, rev = true)
    partialsortperm(ranking, 1:n)
end

struct TradePair
    asset_1::String
    asset_2::String
    original_lt_spread::Float64
    original_st_spread::Float64
    epsilon::Float64
end

function generate_trade_pairs(tickers, st_spreads, lt_spreads, initial_lt_prices, initial_st_prices)
    L_params = calibrate_L.(eachcol(lt_spreads))
    score_L = map(L_params) do (θ, _, σ)
        L_score(θ, σ)
    end
	deltas = map(L_params) do (_, δ, _)
        δ
    end
    params = calibrate.(eachcol(st_spreads), deltas)
    score_S = map(zip(params, deltas)) do ((θ, σ), δ₁)
        score(θ, σ, δ₁)
    end
	pairs = select_pairs(score_L, score_S)
	map(pairs) do idx
        last_spread = st_spreads[end, idx]
        ϵ = quantile(abs.(lt_spreads[2:end, idx] .- lt_spreads[1:end-1, idx]), 0.98)
        pair = index_to_pair(idx)
        TradePair(
            tickers[pair[1]],
            tickers[pair[2]],
            log(initial_lt_prices[pair[1]]) - log(initial_lt_prices[pair[2]]),
            log(initial_st_prices[pair[1]]) - log(initial_st_prices[pair[2]]),
            ϵ 
        )
    end
end
