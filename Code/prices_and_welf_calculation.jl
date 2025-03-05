function Hpindex(fund,L,w,dtradesh;param = param)
    @unpack sigma = param
    a = fund[:,1]
    H = fund[:,2]
    P = (sigma/(sigma-1)).* (w ./ a) .* ((L ./ (sigma .* F .* dtradesh)) .^(1/(1-sigma)))
end
function Hrealw(fund,L,tradesh; param = param)
    @unpack alpha,sigma = param
    a = fund[:,1]
    H = fund[:,2]
    dtradesh = diag(tradesh)
    realwage=((L./(sigma.*F.*dtradesh)).^(alpha./(sigma-1))).*(a.^alpha).*((L./H).^(-(1-alpha)));
    realwage=realwage./(alpha.*((sigma./(sigma-1)).^alpha).*(((1-alpha)./alpha).^(1-alpha)));
end

function Hlandpirce(fund, L, w; param = param)
    @unpack alpha = param
    H = fund[:,2]
    r = ((1-alpha)/alpha) * ((w .* L) ./ H)
end

function Hwelfaregains(ctradesh, tradesh, cL, L; param = param)
    @unpack alpha,sigma = param

    # domestic trade share
    dtradesh = diag(tradesh)
    cdtradesh = diag(ctradesh)

    # welfare gains
    welfgain = ((dtradesh ./ cdtradesh) .^ (alpha / (sigma - 1))) .* ((L ./ cL) .^ (((sigma * (1 - alpha)) - 1) / (sigma - 1)))

    return welfgain
end