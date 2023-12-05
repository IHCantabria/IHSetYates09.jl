
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