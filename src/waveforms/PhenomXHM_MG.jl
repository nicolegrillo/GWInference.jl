"""
PhenomXHM wrappers for +1PN massive-graviton (MG) propagation correction.

`PhenomXHM_TIGER_spinless(1.0)` is the generic TIGER/ppE model:
the last parameter is directly `delta_phi_plus2`.

`PhenomXHM_MG()` is the MG model:
the last parameter is `inv_lambda_g2_km2 = 1 / lambda_g^2` in km^-2,
mapped to the relative +1PN TIGER coefficient `delta_phi_plus2`, using https://arxiv.org/pdf/gr-qc/9709011, Eq. (3.8c, 3.9).
"""

function PolAbs(
    model::PhenomXHM_MG,
    f::AbstractVector,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    inv_lambda_g2_km2,
    Lambda1=0.0,
    Lambda2=0.0,
)
    hp, hc = hphc(model, f, mc, eta, chi1, chi2, dL, iota, inv_lambda_g2_km2)
    return [abs.(hp), abs.(hc)]
end

function Pol(
    model::PhenomXHM_MG,
    f::AbstractVector,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    inv_lambda_g2_km2,
    Lambda1=0.0,
    Lambda2=0.0,
)
    hp, hc = hphc(model, f, mc, eta, chi1, chi2, dL, iota, inv_lambda_g2_km2)
    return [hp, hc]
end

function Phi(
    model::PhenomXHM_MG,
    f,
    mc,
    eta,
    chi1,
    chi2,
    inv_lambda_g2_km2,
    optional_param...;
    fInsJoin_PHI=0.018,
    fcutPar=0.2,
    GMsun_over_c3=uc.GMsun_over_c3,
)
    return f .* 0.0
end

function Ampl(
    model::PhenomXHM_MG,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    inv_lambda_g2_km2;
    fcutPar=0.2,
    fInsJoin_Ampl=0.014,
    GMsun_over_c3=uc.GMsun_over_c3,
    GMsun_over_c2_Gpc=uc.GMsun_over_c2_Gpc,
    GpctoM=uc.uGpc,
)
    hp, hc = hphc(model, f, mc, eta, chi1, chi2, dL, 0.0, inv_lambda_g2_km2)
    return abs.(hp)
end

function _fcut(
    model::PhenomXHM_MG,
    mc,
    eta,
    optional_param...;
    GMsun_over_c3=uc.GMsun_over_c3,
)
    return _fcut(PhenomXHM(), mc, eta; GMsun_over_c3=GMsun_over_c3)
end

function hphc(
    model::PhenomXHM_MG,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    inv_lambda_g2_km2;
    D_Gpc=dL,
    kwargs...,
)
    delta_phi_plus2 = _mg_delta_phi_plus2(mc, eta, D_Gpc, inv_lambda_g2_km2)

    return hphc(
        PhenomXHM_TIGER_spinless(model.PNorder),
        f,
        mc,
        eta,
        chi1,
        chi2,
        dL,
        iota,
        delta_phi_plus2;
        kwargs...,
    )
end

function _mg_beta(mc_Msun, D_Gpc, inv_lambda_g2_km2, GMsun_over_c2=uc.GMsun_over_c2, GpctoM=uc.uGpc)

    "Finds the value of the β parameters for the MG model."

    z = get_z(D_Gpc) # find z to plug into the equation for β

    Mc_km = mc_Msun * GMsun_over_c2 / 1e3 # find the mass in km --> this is the detector mass so when declared it needs a (1+z) term
    D_km = (get_dL_MG(z)[1] / 1e3) * GpctoM / 1e3   # Mpc -> Gpc, then to km

    return pi^2 * D_km * Mc_km * inv_lambda_g2_km2 / (1.0 + z)
end

function _mg_delta_phi_plus2(mc, eta, D_Gpc, inv_lambda_g2_km2)

    "Finds the value of the δϕ^+2 parameters for the MG model."

    beta_MG = _mg_beta(mc, D_Gpc, inv_lambda_g2_km2)
    phi2_GR = (3715.0 / 756.0 + 55.0 * eta / 9.0) * pi^(2.0 / 3.0) # from connection functions

    return -(128.0 / 3.0) * beta_MG * eta^(2.0 / 5.0) * pi^(2/3) / phi2_GR
end

function _mg_lambda_g_bound_km(sigma_inv_lambda_g2_km2)

    return 1.0 / sqrt(sigma_inv_lambda_g2_km2) # find lambda_g in km (if needed) --> could it be useful for after the Fisher? 
end
