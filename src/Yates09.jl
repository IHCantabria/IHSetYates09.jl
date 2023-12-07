function run_Yates09()

    println("Loading libraries...")
    wrkDir = pwd()
    mods = wrkDir*"Modules\\"
    dats = wrkDir*"Data\\"

    # mods = wrkDir*"\\Modules\\"
    # dats = wrkDir*"\\Data\\"

    println("Loading datasets...")

    wavF = NCDataset(dats*"wav.nc")
    ensF = NCDataset(dats*"ens.nc")
    parF = NCDataset(dats*"par.nc")

    configF = NCDataset(dats*"config.nc")

    println("Unpacking datasets...")

    dt = configF["dt"][:][1]
    
    yi = configF["yi"][:][1]

    a = parF["a"][:][1]
    b = parF["b"][:][1]
    cacr = parF["cacr"][:][1]
    cero = parF["cero"][:][1]
    

    
    brk, angBati, depth, Hberm, D50 = configF["brk"][:][1], configF["angBati"][:][1], configF["depth"][:][1], configF["Hberm"][:][1], configF["D50"][:][1]

    if brk == 1
        
        Hs, Tp, θ_w = collect(skipmissing(wavF["Hs"][:])), collect(skipmissing(wavF["Tp"][:])), collect(skipmissing(wavF["Dir"][:]))

        auxAng, auxDepth = similar(Hs), similar(Hs)
        auxAng .= angBati
        auxDepth .= depth

        println("Breaking waves by linear theory...")
        Hb, θ_b, depthb = BreakingPropagation(Hs, Tp, θ_w, auxAng, auxDepth, "spectral")
    else
        Hb, Tp, Hs, depthb = wavF["Hb"][:], wavF["Tp"][:], wavF["Hs"][:], wavF["hb"][:]
    end
    
    close(wavF)
    close(ensF)
    close(configF)
    close(parF)

    println("Datasets closed...")

    Hs = convert(Array{Float64},Hs)
    Tp = convert(Array{Float64},Tp)

    a = convert(Array{Float64},a)
    b = convert(Array{Float64},b)
    cacr = convert(Array{Float64},kacr)
    cero = convert(Array{Float64},kero)

    Yi = convert(Array{Float64},yi)

    ########## START HERE #############

    E = Hb .^ 2

    println("Starting Jaramillo et al. (2021) - Shoreline Rotation Model...")


    Y_t = (E, dt, a, b, cacr, cero, Yi)

    println("\n\n****************Finished****************\n\n")

    return Y_t
    
end

function Yates09(E, dt, a, b, cacr, cero, Yini)
    
    @views Seq = (E .- b) ./ a
    
    Y = similar(E)
    
    Y[1] = Yini

    for i in eachindex(E[2:end])
        if Y[i] < Seq[i+1]
            Y[i+1] = ((Y[i]-Seq[i+1]).*exp(-1. * a *cacr *(E[i+1] ^ 0.5)*dt))+Seq[i+1]
        else
            Y[i+1] = ((Y[i]-Seq[i+1]).*exp(-1. * a *cero *(E[i+1] ^ 0.5)*dt))+Seq[i+1]
        end
    end
    
    return Y, Seq
end

function Yates09Vitousek21Mod(hs, dt, H_s, DY, DT, Yini)
    
    Y = similar(hs)
    
    Y[1] = Yini

    @views Yeq = -DY .*(hs.^2 .- H_s.^2)./H_s.^2
    @views τ = DT .* H_s./hs
    
    for i in eachindex(Y[2:end])
        dYdt = 1 / τ[i+1] * (Yeq[i+1] - Y[i])
        Y[i+1] = Y[i] + dt * dYdt
    end

    return Y, Yeq

end