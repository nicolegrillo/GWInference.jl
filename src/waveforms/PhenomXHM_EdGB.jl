"""
PhenomXHM wrappers for -1PN beyond-GR phase corrections.

`PhenomXHM_TIGER_spinless(-1.0)` is the generic TIGER/ppE model:
the last parameter is directly `delta_phi_minus2`.

`PhenomXHM_EdGB()` is the EdGB model:
the last parameter is `sqrt(alpha_EdGB)` in km, mapped to
`delta_phi_minus2` using arXiv:1905.00870 Eq. (4) and Eq. (10b).
"""

function PolAbs(
    model::PhenomXHM_TIGER_spinless,
    f::AbstractVector,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    delta_phi_minus2,
    Lambda1=0.0,
    Lambda2=0.0,
)
    hp, hc = hphc(model, f, mc, eta, chi1, chi2, dL, iota, delta_phi_minus2)
    return [abs.(hp), abs.(hc)]
end

function Pol(
    model::PhenomXHM_TIGER_spinless,
    f::AbstractVector,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    delta_phi_minus2,
    Lambda1=0.0,
    Lambda2=0.0,
)
    hp, hc = hphc(model, f, mc, eta, chi1, chi2, dL, iota, delta_phi_minus2)
    return [hp, hc]
end

function Phi(
    model::PhenomXHM_TIGER_spinless,
    f,
    mc,
    eta,
    chi1,
    chi2,
    delta_phi_minus2,
    optional_param...;
    fInsJoin_PHI=0.018,
    fcutPar=0.2,
    GMsun_over_c3=uc.GMsun_over_c3,
)
    return f .* 0.0
end

function Ampl(
    model::PhenomXHM_TIGER_spinless,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    delta_phi_minus2;
    fcutPar=0.2,
    fInsJoin_Ampl=0.014,
    GMsun_over_c3=uc.GMsun_over_c3,
    GMsun_over_c2_Gpc=uc.GMsun_over_c2_Gpc,
)
    hp, hc = hphc(model, f, mc, eta, chi1, chi2, dL, 0.0, delta_phi_minus2)
    return abs.(hp)
end

function _fcut(
    model::PhenomXHM_TIGER_spinless,
    mc,
    eta,
    optional_param...;
    GMsun_over_c3=uc.GMsun_over_c3,
)
    return _fcut(PhenomXHM(), mc, eta; GMsun_over_c3=GMsun_over_c3)
end

function hphc(
    model::PhenomXHM_TIGER_spinless,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    delta_phi_minus2;
    kwargs...,
)
    return hphc(
        PhenomXHM(),
        f,
        mc,
        eta,
        chi1,
        chi2,
        dL,
        iota;
        PNorder=model.PNorder,
        o1=delta_phi_minus2,
        kwargs...,
    )
end

function PolAbs(
    model::PhenomXHM_EdGB,
    f::AbstractVector,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    squared_alpha_EdGB_km4,
    Lambda1=0.0,
    Lambda2=0.0,
)
    hp, hc = hphc(model, f, mc, eta, chi1, chi2, dL, iota, squared_alpha_EdGB_km4)
    return [abs.(hp), abs.(hc)]
end

function Pol(
    model::PhenomXHM_EdGB,
    f::AbstractVector,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    squared_alpha_EdGB_km4,
    Lambda1=0.0,
    Lambda2=0.0,
)
    hp, hc = hphc(model, f, mc, eta, chi1, chi2, dL, iota, squared_alpha_EdGB_km4)
    return [hp, hc]
end

function Phi(
    model::PhenomXHM_EdGB,
    f,
    mc,
    eta,
    chi1,
    chi2,
    squared_alpha_EdGB_km4,
    optional_param...;
    fInsJoin_PHI=0.018,
    fcutPar=0.2,
    GMsun_over_c3=uc.GMsun_over_c3,
)
    return f .* 0.0
end

function Ampl(
    model::PhenomXHM_EdGB,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    squared_alpha_EdGB_km4;
    fcutPar=0.2,
    fInsJoin_Ampl=0.014,
    GMsun_over_c3=uc.GMsun_over_c3,
    GMsun_over_c2_Gpc=uc.GMsun_over_c2_Gpc,
)
    hp, hc = hphc(model, f, mc, eta, chi1, chi2, dL, 0.0, squared_alpha_EdGB_km4)
    return abs.(hp)
end

function _fcut(
    model::PhenomXHM_EdGB,
    mc,
    eta,
    optional_param...;
    GMsun_over_c3=uc.GMsun_over_c3,
)
    return _fcut(PhenomXHM(), mc, eta; GMsun_over_c3=GMsun_over_c3)
end

function hphc(
    model::PhenomXHM_EdGB,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    squared_alpha_EdGB_km4;
    kwargs...,
)
    delta_phi_minus2 = _edgb_delta_phi_minus2(mc, eta, chi1, chi2, squared_alpha_EdGB_km4)

    return hphc(
        PhenomXHM_TIGER_spinless(model.PNorder),
        f,
        mc,
        eta,
        chi1,
        chi2,
        dL,
        iota,
        delta_phi_minus2;
        kwargs...,
    )
end

function _edgb_scalar_charge_factor(chi)
    if abs(chi) < 1e-4
        return 1.0 - chi^2 / 4.0
    end

    return 2.0 * (sqrt(1.0 - chi^2) - 1.0 + chi^2) / chi^2
end

function _edgb_delta_phi_minus2(mc, eta, chi1, chi2, squared_alpha_EdGB_km4, GMsun_over_c2=uc.GMsun_over_c2)

    GMsun_over_c2_km = GMsun_over_c2 / 1e3 # convert to km
    total_mass_km = mc / eta^(3.0 / 5.0) * GMsun_over_c2_km
    root = ifelse(eta < 0.25, sqrt(1.0 - 4.0 * eta), 0.0)
    # find two masses both in km
    m1 = 0.5 * total_mass_km * (1.0 + root)
    m2 = 0.5 * total_mass_km * (1.0 - root)

    s1 = _edgb_scalar_charge_factor(chi1)
    s2 = _edgb_scalar_charge_factor(chi2)

    zeta_EdGB = 16.0 * pi * squared_alpha_EdGB_km4 / total_mass_km^4

    beta_EdGB =
        -5.0 / 7168.0 *
        zeta_EdGB *
        (m1^2 * s2 - m2^2 * s1)^2 /
        (total_mass_km^4 * eta^(18.0 / 5.0))

    return (128.0 / 3.0) * beta_EdGB * eta^(-2.0 / 5.0) 
end

function _edgb_small_coupling_ratio(mc, eta, sqrt_alpha_EdGB_km, GMsun_over_c2=uc.GMsun_over_c2)
    
    GMsun_over_c2_km = GMsun_over_c2 / 1e3 # convert to km
    total_mass_km = mc / eta^(3.0 / 5.0) * GMsun_over_c2_km
    root = ifelse(eta < 0.25, sqrt(1.0 - 4.0 * eta), 0.0)
    m1 = 0.5 * total_mass_km * (1.0 + root)
    m2 = 0.5 * total_mass_km * (1.0 - root)

    return sqrt_alpha_EdGB_km / min(m1, m2)
end
