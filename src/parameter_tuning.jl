# Estimate Y-process

function variance(θ, δ, A, B, C, D, E, F)
    a = exp(-θ*δ)
    √(2/(77*30)*(θ/(1-a^2)*(A*a^2+B*a+C) - θ/(1-a^156)*(D*a^156+E*a^78+F)))
end

function log_likelihood(Y, L, δ, A, B, C, D, E, F, θ)
    a = exp(-θ*δ)
    σ = variance(θ, δ, A, B, C, D, E, F)
    -77 * 30log(σ)-77*30/2*log((1-a^2)/θ)-77*30/2*log(π)+30/2*log((1-a^156)/(1-a^2))-77*30/2
end

function calibrate(Y, δ₁)
    try
        L = zeros(60)
        for i in 1:30
            L[2i-1] = Y[79*(i-1)+1]
            L[2i] = Y[79*i]
        end
        L = ShiftedArray(vcat(0, L), -1)
        δ = δ₁/78
        A = sum(sum(Y[79*(i-1)+j]^2 for j in 1:78) + 78*(L[2i-2] + L[2i-1])^2/4 - (L[2i-2] + L[2i-1])*sum(Y[79*(i-1)+j] for j in 1:78) for i in 1:30)
        B = sum(-156/4*(L[2i-2]+L[2i-1])^2 + (L[2i-2] + L[2i-1]) * sum(Y[79*(i-1)+j+1] for j in 1:78) + (L[2i-2] + L[2i-1]) * sum(Y[79*(i-1)+j] for j in 1:78) - 2 * sum(Y[79*(i-1)+j+1]*Y[79*(i-1)+j] for j in 1:78) for i in 1:30)
        C = sum(sum(Y[79*(i-1)+j+1]^2 for j in 1:78) + 78/4* (L[2i-2] + L[2i-1])^2 - (L[2i-2] + L[2i-1])* sum(Y[79*(i-1)+j+1] for j in 1:78) for i in 1:30)
        D = sum(L[2i-1]^2 / 4 + L[2i-2]^2 / 4 - L[2i-1]*L[2i-2]/2 for i in 1:30)
        E = sum(L[2i-1]^2/2 - L[2i-2]^2/2 - L[2i]*L[2i-1] + L[2i]*L[2i-2] for i in 1:30)
        F = sum(L[2i]^2 + L[2i-1]^2 / 4 + L[2i-2]^2 / 4 - L[2i]*L[2i-1] - L[2i]*L[2i-2] + L[2i-1]*L[2i-2]/2 for i in 1:30)
        opt = θ -> -log_likelihood(Y, L, δ, A, B, C, D, E, F, only(θ))
        init = [0.01]
        solve = optimize(opt, init, BFGS())
        θ = only(minimizer(solve))
        σ = variance(θ, δ, A, B, C, D, E, F)
        (θ, σ)
    catch e
        (-Inf, -Inf)
    end
end

# Estimate L-process

@views function intraday_returns(L)
    L[2:2:end] .- L[1:2:end]
end

@views function overnight_returns(L)
    L[3:2:end] .- L[2:2:end-2]
end

function initialize_deltas(L)
    intraday = intraday_returns(L)
    overnight = overnight_returns(L)
    ratio = var(intraday) / var(overnight)
    δ₂ = 1.0/(252.0 * (ratio + 1.0))
    δ₁ = 1.0/252.0 - δ₂
    (δ₁, δ₂)
end

@views function log_likelihood(L, δ₁, δ₂, θ, σ₁, σ₂)
    (σ₁ < 0 || σ₂ < 0) && return -Inf
    length(L) % 2 == 0 || error("not even")
    N = length(L)÷2
    -N / 2 * log(2π) - N * log(σ₁) - 1/(2σ₁^2)*sum((L[2i]-L[2i-1]*exp(-θ*δ₁))^2 for i in 1:N) - (N-1)/2 * log(2π) - (N-1)*log(σ₂) - 1/(2σ₂^2)*sum((L[2i+1]-L[2i]*exp(-θ*δ₂))^2 for i in 1:N-1)
end

@views function refine_deltas(L, θ, δ₁::Float64)
    f₁ = δ₁ -> var(L[2:2:end] .- L[1:2:end]*exp(-θ*δ₁)) / var(L[3:2:end] .- L[2:2:end-2]*exp(-θ*(1/252 - δ₁))) - (1-exp(-2θ*δ₁))/(1-exp(-2θ*(1/252-δ₁)))
    δ₁ = find_zero(f₁, δ₁)
    δ₂ = 1/252 - δ₁
    δ₁, δ₂
end

function calibrate_L(L)
    intraday = intraday_returns(L)
    overnight = overnight_returns(L)
    tol::Float64 = Inf
    δ₁, δ₂ = initialize_deltas(L)
    θ, σ₁, σ₂ = 0.0, 0.0, 0.0
    iterations = 0
    while tol > 1e-6
        if iterations > 50
            #@warn "Max iterations reached. Tolerance: $tol"
            break
        end
        opt = params -> -log_likelihood(L, δ₁, δ₂, params...)
        init = [0.0, 0.1, 0.1]
        lower = [-Inf, 0.0, 0.0]
        upper = [Inf, Inf, Inf]
        solve = optimize(opt, lower, upper, init, Fminbox(BFGS()))
        θ, σ₁, σ₂ = minimizer(solve)
        δ₁_temp, δ₂_temp = refine_deltas(L, θ, δ₁)
        tol = abs(δ₁ - δ₁_temp) + abs(δ₂ - δ₂_temp)
        δ₁, δ₂ = δ₁_temp, δ₂_temp
        iterations += 1
    end
    σ = σ₁ / √((1-exp(-2θ*δ₁))/(2θ))
    θ, δ₁, σ
end

