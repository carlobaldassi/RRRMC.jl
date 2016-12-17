module Newton

export newton, NewtonParameters

function ∇!(∂f::Matrix, f::Function, x0, δ, f0, x1)
    n = length(x0)
    copy!(x1, x0)
    for i = 1:n
        x1[i] += δ
        ∂f[:,i] = (f(x1) - f0) / δ
        x1[i] = x0[i]
    end
end

∇(f::Function, x0::Real, δ::Real, f0::Real) = (f(x0 + δ) - f0) / δ

type NewtonParameters
    δ::Float64
    ϵ::Float64
    verb::Int
    maxiters::Int
end

function newton(f::Function, x₀::Real, pars::NewtonParameters)
    η = 1.0
    ∂f = 0.0
    x = x₀
    x1 = 0.0

    f0 = f(x)
    @assert isa(f0, Real)
    normf0 = abs(f0)
    it = 0
    while normf0 ≥ pars.ϵ
        it > pars.maxiters && return (false, x, it, normf0)
        it += 1
        if pars.verb > 1
            println("(𝔫) it=$it")
            println("(𝔫)   x=$x")
            println("(𝔫)   f0=$f0")
            println("(𝔫)   norm=$(abs(f0))")
            println("(𝔫)   η=$η")
        end
        δ = pars.δ
        while true
            try
                ∂f = ∇(f, x, δ, f0)
                break
            catch
                δ /= 2
            end
            if δ < 1e-15
                normf0 ≥ pars.ϵ && warn("newton:  δ=$δ")
                return (false, x, it, normf0)
            end
        end
        Δx = -f0 / ∂f
        pars.verb > 1 && println("(𝔫)  Δx=$Δx")
        while true
            x1 = x + Δx * η
            local new_f0, new_normf0
            try
                new_f0 = f(x1)
                new_normf0 = abs(new_f0)
            catch
                new_normf0 = Inf
            end
            if new_normf0 < normf0
                η = min(1.0, η * 1.1)
                f0 = new_f0
                normf0 = new_normf0
                x = x1
                break
            end
            η /= 2
            η < 1e-15 && return (false, x, it, normf0)
        end
    end
    return true, x, it, normf0
end

function newton(f::Function, x₀, pars::NewtonParameters)
    η = 1.0
    n = length(x₀)
    ∂f = Array(Float64, n, n)
    x = Float64[x₀[i] for i = 1:n]
    x1 = Array(Float64, n)

    f0 = f(x)
    @assert length(f0) == n
    @assert isa(f0, Union{Real,Vector})
    normf0 = vecnorm(f0)
    it = 0
    while normf0 ≥ pars.ϵ
        it > pars.maxiters && return (false, x, it, normf0)
        it += 1
        if pars.verb > 1
            println("(𝔫) it=$it")
            println("(𝔫)   x=$x")
            println("(𝔫)   f0=$f0")
            println("(𝔫)   norm=$(vecnorm(f0))")
            println("(𝔫)   η=$η")
        end
        δ = pars.δ
        while true
            try
                ∇!(∂f, f, x, δ, f0, x1)
                break
            catch
                δ /= 2
            end
            if δ < 1e-15
                normf0 ≥ pars.ϵ && warn("newton:  δ=$δ")
                return (false, x, it, normf0)
            end
        end
        if isa(f0, Vector)
            Δx = -∂f \ f0
        else
            Δx = -f0 / ∂f[1,1]
        end
        pars.verb > 1 && println("(𝔫)  Δx=$Δx")
        while true
            for i = 1:n
                x1[i] = x[i] + Δx[i] * η
            end
            local new_f0, new_normf0
            try
                new_f0 = f(x1)
                new_normf0 = vecnorm(new_f0)
            catch
                new_normf0 = Inf
            end
            if new_normf0 < normf0
                η = min(1.0, η * 1.1)
                if isa(f0, Vector)
                    copy!(f0, new_f0)
                else
                    f0 = new_f0
                end
                normf0 = new_normf0
                copy!(x, x1)
                break
            end
            η /= 2
            η < 1e-15 && return (false, x, it, normf0)
        end
    end
    return true, x, it, normf0
end

end # module
