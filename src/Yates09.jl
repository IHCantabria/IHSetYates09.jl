"""
#######################################
# Yates et al. (2004) model           #
# E: wave energy                      #
# dt: time step                       #
# a: model's parameter                #
# b: model's parameter                #
# kero: erosion coefficient           #
# kacr: accretion coefficient         #
# Yi: initial position                #
#######################################
#######################################
# Vitousek et al. (2021) model        #
# Hs: wave significant height         #
# dt: time step                       #
# H_s: model's parameter              #
# DY: model's parameter               #
# DT: model's parameter               #
# Yi: initial position                #
#######################################
"""

function run_Yates09()

    println("Loading libraries...")
    wrkDir = pwd()

    println("Loading datasets...")

    wavF = dats*"wav.nc"
    parF = dats*"par.nc"
    slF = dats*"sl.nc"
    configF = dats*"config.nc"

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
    
    Seq = (E .- b) ./ a
    
    Y = similar(E)
    
    Y[1] = Yini

    for i in eachindex(E[2:end])
        if Y[i] < Seq[i+1]
            Y[i+1] = ((Y[i]-Seq[i+1]).*exp(-1. * a *cacr *(E[i+1] ^ 0.5)*dt))+Seq[i+1]
        else
            Y[i+1] = ((Y[i]-Seq[i+1]).*exp(-1. * a *cero *(E[i+1] ^ 0.5)*dt))+Seq[i+1]
        end
    end
    
    return Y
end

function Yates09Vitousek21Mod(hs, dt, H_s, DY, DT, Yini)
    
    Y = zeros(length(hs))
    
    Y[1] = Yini

    @views Yeq = -DY .*(hs.^2 .- H_s.^2)./H_s.^2
    @views τ = DT .* H_s./hs
    
    for i in eachindex(Y[2:end])
        dYdt = 1 / τ[i+1] * (Yeq[i+1] - Y[i])
        Y[i+1] = Y[i] + dt * dYdt
    end

    return Y, Yeq

end

function cal_Yates09()

    println("Loading libraries...")
    wrkDir = pwd()
    dats = wrkDir*"/data/"

    # mods = wrkDir*"\\Modules\\"
    # dats = wrkDir*"\\Data\\"

    println("Reading datasets...")

    wavF = dats*"wav.nc"
    ensF = dats*"ens.nc"
    parF = dats*"par.nc"

    configF = dats*"config.nc"

    dt = ncread(configF, "dt")[1]
    calPar = ncread(configF, "calPar")[1]
    
    brk, angBati, depth, D50 = ncread(configF, "brk")[1], ncread(configF, "angBati")[1], ncread(configF, "depth")[1], ncread(configF, "D50")[1]

    MetObj = ncread(configF, "MetObj")[1]

    if brk == 1
        
        Hs, Tp, θ_w = ncread(wavF, "Hs"), ncread(wavF, "Tp"), ncread(wavF, "Dir")

        auxAng, auxDepth = similar(Hs), similar(Hs)
        auxAng .= angBati
        auxDepth .= depth

        println("Breaking waves by linear theory...")
        Hb, θ_b, depthb = BreakingPropagation(Hs, Tp, θ_w, auxAng, auxDepth, "spectral")
    else
        Hb, Tp, Hs, depthb = ncread(wavF, "Hb"), ncread(wavF, "Tp"), ncread(wavF, "Hs"), ncread(wavF, "hb")
    end

    YY, MM, DD, HH = ncread(wavF, "Y"), ncread(wavF, "M"), ncread(wavF, "D"), ncread(wavF, "h")

    YYo, MMo, DDo, HHo = ncread(ensF, "Y"), ncread(ensF, "M"), ncread(ensF, "D"), ncread(ensF, "h")
    
    Y_obs = ncread(ensF, "Obs")

    t_obs = DateTime.(YYo, MMo, DDo, HHo)
    t_wav = DateTime.(YY, MM, DD, HH)

    ii =  t_obs .<= t_wav[end] .&& t_obs .>= t_wav[1]

    t_obs, Y_obs = t_obs[ii], Y_obs[ii]

    ii =  t_wav .<= t_obs[end] .&& t_wav .>= t_obs[1]

    t_wav, hb, tp, hs, depthb = t_wav[ii], Hb[ii], Tp[ii], Hs[ii], depthb[ii]

    idx_obs = zeros(length(t_obs))

    for i in eachindex(t_obs)
        idx_obs[i] = argmin(abs.(t_wav .- t_obs[i]))
    end

    idx_obs = convert(Array{Int64},idx_obs)

    ########## START HERE #############

    E = Hb .^ 2

    println("Starting Yates et al. (2009) - Shoreline Evolution Model...")

    Ymdr = Dict()
    aRP = Dict()
    aRMSE = Dict()
    aMSS = Dict()
    popr = Dict()
    objr = Dict()

    if calPar == 4
        println("calPar =" * string(calPar) * " - 4 parameters")
        function Calibra_4r(Χ)
            println("calPar =" * string(calPar) * " 214")
            Ymd = Yates09(E, dt, -exp(Χ[1]), exp(Χ[2]), -exp(Χ[3]), -exp(Χ[4]), Y_obs[1])
            YYsl = Ymd[idx_obs]
            if MetObj == "Pearson"
                return 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl)))
            elseif MetObj == "RMSE"
                return abs(sqrt(mean((YYsl .- Y_obs).^2))/5)
            elseif MetObj == "MSS"
                return sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2)
            elseif MetObj == "BSS"
                return (mean((YYsl .- Y_obs).^2) - mean((YYref .- Y_obs).^2))/mean((YYref .- Y_obs).^2)
            elseif MetObj == "Double"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5))
            elseif MetObj == "Triple"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double2"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double3"
                return (abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            end
        end

        
        boundsr = [(log(1e-3), log(5e-1)),
                    (log(1e-1), log(1e+2)),
                    (log(1e-5), log(1e-1)),
                    (log(1e-5), log(1e-1))]

        if MetObj == "Double" || MetObj == "Double2" || MetObj == "Double3"
            println("calPar =" * string(calPar) * " 243")
            resr = bboptimize(Calibra_4r; 
                            # Method = :simultaneous_perturbation_stochastic_approximation,
                            SearchRange = boundsr,
                            NumDimensions = 4,
                            PopulationSize = 500,
                            MaxSteps = 5000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 10000,
                            Method=:borg_moea)
        elseif MetObj == "Triple"
            resr = bboptimize(Calibra_4r; 
                            # Method = :simultaneous_perturbation_stochastic_approximation,
                            SearchRange = boundsr,
                            NumDimensions = 4,
                            PopulationSize = 500,
                            MaxSteps = 5000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{3}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 10000,
                            Method=:borg_moea)
        else
            resr = bboptimize(Calibra_4r; 
                            Method = :adaptive_de_rand_1_bin,
                            SearchRange = boundsr,
                            NumDimensions = 4,
                            PopulationSize = 500,
                            MaxSteps = 5000,
                            FitnessTolerance = 1e-6,
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 10000)
        end

        objr = best_fitness(resr)
        popr = best_candidate(resr)

        Ymdr = Yates09(E, dt, -exp(popr[1]), exp(popr[2]), -exp(popr[3]), -exp(popr[4]), Y_obs[1])

        Ysl = Ymdr[idx_obs]
        aRP = sum((Ysl.-mean(Ysl)).*(Y_obs .- mean(Y_obs)))/(std(Ysl)*std(Y_obs)*length(Ysl))
        aRMSE = sqrt(mean((Ysl .- Y_obs).^2))
        aMSS = 1 - sum((Ysl .- Y_obs).^2)/length(Ysl)/(var(Ysl)+var(Y_obs)+(mean(Ysl)-mean(Y_obs))^2)


        println("\n\n****************Finished****************\n\n")

        year_atts = Dict("long_name" => "Year")
        month_atts = Dict("long_name" => "Month")
        day_atts = Dict("long_name" => "Day")
        hour_atts = Dict("long_name" => "Hour")
        println("Writing output...")

        output = wrkDir*"/results/Shoreline_Y09.nc"
        nccreate(output, "year",
                    "dim", length(YY),
                    atts = year_atts)
        ncwrite(YY, output, "year")
        nccreate(output, "month",
                    "dim", length(MM),
                    atts = month_atts)
        ncwrite(MM, output, "month")
        nccreate(output, "day",
                    "dim", length(DD),
                    atts = day_atts)
        ncwrite(DD, output, "day")
        nccreate(output, "hour",
                    "dim", length(HH),
                    atts = hour_atts)
        ncwrite(HH, output, "hour")  

        Y_atts = Dict("units" => "m",
            "long_name" => "Shoreline position",
            "standard_name" => "Y")
        Yi_atts = Dict("units" => "m",
            "long_name" => "Initial shoreline position",
            "standard_name" => "Yi")
        kacr_atts = Dict("units" => "-",
            "long_name" => "Accretion coefficient",
            "standard_name" => "cacr")
        kero_atts = Dict("units" => "-",
            "long_name" => "Erosion coefficient",
            "standard_name" => "cero")
        a_atts = Dict("units" => "-",
            "long_name" => "a parameter",
            "standard_name" => "a")
        b_atts = Dict("units" => "-",
            "long_name" => "b parameter",
            "standard_name" => "b")
        RP_atts = Dict("units" => "-",
            "long_name" => "Pearson correlation coefficient",
            "standard_name" => "RP")
        RMSE_atts = Dict("units" => "m",
            "long_name" => "Root mean square error",
            "standard_name" => "RMSE")
        MSS_atts = Dict("units" => "-",
            "long_name" => "Mielke Skill Score",
            "standard_name" => "MSS")


        nccreate(output, "Y",
                    "dim", length(Ymdr),
                    atts = Y_atts)
        ncwrite(Ymdr, output, "Y")
        nccreate(output, "Yi",
                    "len", 1,
                    atts = Yi_atts)
        ncwrite([Y_obs[1]], output, "Yi")
        nccreate(output, "a",
                    "len", 1,
                    atts = a_atts)
        ncwrite([exp(popr[1])], output, "a")
        nccreate(output, "b",
                    "len", 1,
                    atts = b_atts)
        ncwrite([exp(popr[2])], output, "b")
        nccreate(output, "cacr",
                    "len", 1,
                    atts = kacr_atts)
        ncwrite([exp(popr[3])], output, "cacr")
        nccreate(output, "cero",
                    "len", 1,
                    atts = kero_atts)
        ncwrite([exp(popr[4])], output, "cero")
        nccreate(output, "RP",
                    "len", 1,
                    atts = RP_atts)
        ncwrite([aRP], output, "RP")
        nccreate(output, "RMS",
                    "len", 1,
                    atts = RMSE_atts)
        ncwrite([aRMSE], output, "RMSE_flagP="*string(i))
        nccreate(output, "MSS_flagP="*string(i),
                    "len", 1,
                    atts = MSS_atts)
        ncwrite([aMSS], output, "MSS_flagP="*string(i))

    elseif calPar == 5
        println("calPar =" * string(calPar) * " - 388")
        function Calibra_5r(Χ)
            println("calPar =" * string(calPar) * " 390")
            Ymd = Yates09(E, dt, -exp(Χ[1]), exp(Χ[2]), -exp(Χ[3]), -exp(Χ[4]), Χ[5])
            YYsl = Ymd[idx_obs]
            if MetObj == "Pearson"
                return 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl)))
            elseif MetObj == "RMSE"
                return abs(sqrt(mean((YYsl .- Y_obs).^2))/5)
            elseif MetObj == "MSS"
                return sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2)
            elseif MetObj == "BSS"
                return (mean((YYsl .- Y_obs).^2) - mean((YYref .- Y_obs).^2))/mean((YYref .- Y_obs).^2)
            elseif MetObj == "Double"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5))
            elseif MetObj == "Triple"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double2"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            elseif MetObj == "Double3"
                return (abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            end
        end

        
        boundsr = [(log(1e-3), log(5e-1)),
                    (log(1e-1), log(1e+2)),
                    (log(1e-5), log(1e-1)),
                    (log(1e-5), log(1e-1)),
                    (0.25*minimum(Y_obs), 2*maximum(Y_obs))]

        if MetObj == "Double" || MetObj == "Double2" || MetObj == "Double3"
            println("calPar =" * string(calPar) * " 420")
            resr = bboptimize(Calibra_5r; 
                            # Method = :simultaneous_perturbation_stochastic_approximation,
                            SearchRange = boundsr,
                            NumDimensions = 5,
                            PopulationSize = 500,
                            MaxSteps = 5000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 10000,
                            Method=:borg_moea)
        elseif MetObj == "Triple"
            resr = bboptimize(Calibra_5r; 
                            # Method = :simultaneous_perturbation_stochastic_approximation,
                            SearchRange = boundsr,
                            NumDimensions = 5,
                            PopulationSize = 500,
                            MaxSteps = 5000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{3}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 10000,
                            Method=:borg_moea)
        else
            resr = bboptimize(Calibra_5r; 
                            Method = :adaptive_de_rand_1_bin,
                            SearchRange = boundsr,
                            NumDimensions = 5,
                            PopulationSize = 500,
                            MaxSteps = 5000,
                            FitnessTolerance = 1e-6,
                            TraceMode=:compact,
                            ϵ=0.1,
                            τ = 0.05,
                            MaxStepsWithoutEpsProgress = 10000)
        end

        objr = best_fitness(resr)
        popr = best_candidate(resr)

        Ymdr = Yates09(E, dt, -exp(popr[1]), exp(popr[2]), -exp(popr[3]), -exp(popr[4]), popr[5])

        Ysl = Ymdr[idx_obs]
        aRP = sum((Ysl.-mean(Ysl)).*(Y_obs .- mean(Y_obs)))/(std(Ysl)*std(Y_obs)*length(Ysl))
        aRMSE = sqrt(mean((Ysl .- Y_obs).^2))
        aMSS = 1 - sum((Ysl .- Y_obs).^2)/length(Ysl)/(var(Ysl)+var(Y_obs)+(mean(Ysl)-mean(Y_obs))^2)

        println("\n\n****************Writing output****************\n\n")

        year_atts = Dict("long_name" => "Year")
        month_atts = Dict("long_name" => "Month")
        day_atts = Dict("long_name" => "Day")
        hour_atts = Dict("long_name" => "Hour")
        println("Writing output...")

        output = wrkDir*"/results/Shoreline_Y09.nc"
        nccreate(output, "year",
                    "dim", length(YY),
                    atts = year_atts)
        ncwrite(YY, output, "year")
        nccreate(output, "month",
                    "dim", length(MM),
                    atts = month_atts)
        ncwrite(MM, output, "month")
        nccreate(output, "day",
                    "dim", length(DD),
                    atts = day_atts)
        ncwrite(DD, output, "day")
        nccreate(output, "hour",
                    "dim", length(HH),
                    atts = hour_atts)
        ncwrite(HH, output, "hour")  

        Y_atts = Dict("units" => "m",
            "long_name" => "Shoreline position",
            "standard_name" => "Y")
        Yi_atts = Dict("units" => "m",
            "long_name" => "Initial shoreline position",
            "standard_name" => "Yi")
        kacr_atts = Dict("units" => "-",
            "long_name" => "Accretion coefficient",
            "standard_name" => "cacr")
        kero_atts = Dict("units" => "-",
            "long_name" => "Erosion coefficient",
            "standard_name" => "cero")
        a_atts = Dict("units" => "-",
            "long_name" => "a parameter",
            "standard_name" => "a")
        b_atts = Dict("units" => "-",
            "long_name" => "b parameter",
            "standard_name" => "b")
        RP_atts = Dict("units" => "-",
            "long_name" => "Pearson correlation coefficient",
            "standard_name" => "RP")
        RMSE_atts = Dict("units" => "m",
            "long_name" => "Root mean square error",
            "standard_name" => "RMSE")
        MSS_atts = Dict("units" => "-",
            "long_name" => "Mielke Skill Score",
            "standard_name" => "MSS")

        nccreate(output, "Y",
                    "dim", length(Ymdr),
                    atts = Y_atts)
        ncwrite(Ymdr, output, "Y")
        nccreate(output, "Yi",
                    "len", 1,
                    atts = Yi_atts)
        ncwrite([popr[5]], output, "Yi")
        nccreate(output, "a",
                    "len", 1,
                    atts = a_atts)
        ncwrite([exp(popr[1])], output, "a")
        nccreate(output, "b",
                    "len", 1,
                    atts = b_atts)
        ncwrite([exp(popr[2])], output, "b")
        nccreate(output, "cacr",
                    "len", 1,
                    atts = kacr_atts)
        ncwrite([exp(popr[3])], output, "cacr")
        nccreate(output, "cero",
                    "len", 1,
                    atts = kero_atts)
        ncwrite([exp(popr[4])], output, "cero")
        nccreate(output, "RP",
                    "len", 1,
                    atts = RP_atts)
        ncwrite([aRP], output, "RP")
        nccreate(output, "RMS",
                    "len", 1,
                    atts = RMSE_atts)
        ncwrite([aRMSE], output, "RMSE_flagP="*string(i))
        nccreate(output, "MSS_flagP="*string(i),
                    "len", 1,
                    atts = MSS_atts)
        ncwrite([aMSS], output, "MSS_flagP="*string(i))

    end

    println("\n\n****************Finished****************\n\n")

end