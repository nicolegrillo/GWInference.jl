"""
PhenomXHM wrappers for 2PN beyond-GR phase corrections.

`PhenomXHM_dCS()` is the dCS model:
the last parameter is `alpha_dCS^2` in km^4, mapped to
`delta_phi_plus4` using arXiv:1905.00870 Eq. (3).
"""

function PolAbs(
    model::PhenomXHM_dCS,
    f::AbstractVector,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    squared_alpha_dCS_km4,
    Lambda1=0.0,
    Lambda2=0.0,
)
    hp, hc = hphc(model, f, mc, eta, chi1, chi2, dL, iota, squared_alpha_dCS_km4)
    return [abs.(hp), abs.(hc)]
end

function Pol(
    model::PhenomXHM_dCS,
    f::AbstractVector,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    squared_alpha_dCS_km4,
    Lambda1=0.0,
    Lambda2=0.0,
)
    hp, hc = hphc(model, f, mc, eta, chi1, chi2, dL, iota, squared_alpha_dCS_km4)
    return [hp, hc]
end

function Phi(
    model::PhenomXHM_dCS,
    f,
    mc,
    eta,
    chi1,
    chi2,
    squared_alpha_dCS_km4,
    optional_param...;
    fInsJoin_PHI=0.018,
    fcutPar=0.2,
    GMsun_over_c3=uc.GMsun_over_c3,
)
    return f .* 0.0
end

function Ampl(
    model::PhenomXHM_dCS,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    squared_alpha_dCS_km4;
    fcutPar=0.2,
    fInsJoin_Ampl=0.014,
    GMsun_over_c3=uc.GMsun_over_c3,
    GMsun_over_c2_Gpc=uc.GMsun_over_c2_Gpc,
)
    hp, hc = hphc(model, f, mc, eta, chi1, chi2, dL, 0.0, squared_alpha_dCS_km4)
    return abs.(hp)
end

function _fcut(
    model::PhenomXHM_dCS,
    mc,
    eta,
    optional_param...;
    GMsun_over_c3=uc.GMsun_over_c3,
)
    return _fcut(PhenomXHM(), mc, eta; GMsun_over_c3=GMsun_over_c3)
end

function hphc(
    model::PhenomXHM_dCS,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    squared_alpha_dCS_km4;
    kwargs...,
)
    delta_phi_plus4 = _dcs_delta_phi_plus4(mc, eta, chi1, chi2, squared_alpha_dCS_km4)

    return hphc(
        PhenomXHM_TIGER_spinless(model.PNorder),
        f,
        mc,
        eta,
        chi1,
        chi2,
        dL,
        iota,
        delta_phi_plus4;
        kwargs...,
    )
end

function _dcs_scalar_charge_factor(χ)
    if abs(χ) < 1e-4
        return 5.0 * χ / 8.0
    end

    return (2 + 2 * χ^4 - 2 * sqrt(1 - χ^2) - χ^2 * (3 - 2 * sqrt(1 - χ^2))) / (2 * χ^3)
end

function _dcs_delta_phi_plus4(mc, eta, chi1, chi2, squared_alpha_dCS_km4, GMsun_over_c2=uc.GMsun_over_c2)

    GMsun_over_c2_km = GMsun_over_c2 / 1e3 # convert to km
    total_mass_km = mc / eta^(3.0 / 5.0) * GMsun_over_c2_km
    root = ifelse(eta < 0.25, sqrt(1.0 - 4.0 * eta), 0.0)
    # find two masses both in km
    m1 = 0.5 * total_mass_km * (1.0 + root)
    m2 = 0.5 * total_mass_km * (1.0 - root)

    s1 = _dcs_scalar_charge_factor(chi1)
    s2 = _dcs_scalar_charge_factor(chi2)

    ζ_dCS = 16.0 * pi * squared_alpha_dCS_km4 / total_mass_km^4

    η = eta
    eta2 = eta^2
    delta = sqrt(max(0.0, 1.0 - 4.0 * eta))
    chi12 = chi1^2
    chi22 = chi2^2
    chi1dotchi2 = chi1 * chi2
    phi4_GR = ((15293365. /508032. + (27145. *eta)/504. + (3085*eta2)/72.)*(pi^(4. /3.)) + ((-5. *(81. *chi12*(1. + delta - 2. *eta) + 316. *chi1dotchi2*eta - 81. *chi22*(-1. + delta + 2. *eta)))/16.)*(pi^(4. /3.)))

    β_dCS =
        -5.0 / 8192.0 * ζ_dCS / η^(14. / 5.) * 
        (m1 * s2 - m2 * s1)^2 / total_mass_km^2 + 
        15075.0 / 114688.0 * ζ_dCS / η^(14.0 / 5.0) * 1 / total_mass_km^2 * 
        (m2^2 * chi1^2 - (350.0 / 201.0) * m1 * m2 * chi1 * chi2 + m1^2 * chi2^2)


    return (128 / 3) * β_dCS * η^(4/5) / phi4_GR # not so sure about conversion
end

function _dcs_small_coupling_ratio(mc, eta, sqrt_alpha_dCS_km, GMsun_over_c2=uc.GMsun_over_c2)

    GMsun_over_c2_km = GMsun_over_c2 / 1e3 # convert to km
    total_mass_km = mc / eta^(3.0 / 5.0) * GMsun_over_c2_km
    root = ifelse(eta < 0.25, sqrt(1.0 - 4.0 * eta), 0.0)
    m1 = 0.5 * total_mass_km * (1.0 + root)
    m2 = 0.5 * total_mass_km * (1.0 - root)

    return sqrt_alpha_dCS_km / min(m1, m2)
end
