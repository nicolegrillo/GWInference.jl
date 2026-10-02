"""
PhenomXHM wrapper for the eccentricity phase correction, to O(e0^8) in the
leading (Newtonian point-particle, dominant ℓ=2 harmonic) eccentric expansion
of Yunes, Arun, Berti & Will, arXiv:0906.0313, Eq. (4.28).

Unlike the EdGB/dCS/MG corrections, eccentricity is NOT routed through the
generic `PNorder`/`o1` ppE mechanism: Eq. (4.28) needs 4 simultaneous new PN
order terms (-19/6, -19/3, -19/2, -38/3PN, see `_ecc_phase_coeffs` in
ConnectionFunctionsXAS.jl), not a single new/modified coefficient, so it
can't be selected via a single PNorder switch. Instead `e0`/`f0` are their
own keyword arguments, independent of PNorder/o1, threaded directly through
`hphc(model::PhenomXHM, ...)`. This also means eccentricity composes freely
with any of EdGB/dCS/MG: e.g. `hphc(PhenomXHM_TIGER_spinless(-1.0), ...,
delta_phi_minus2; e0=e0, f0=f0)` would give EdGB + eccentricity together.

`PhenomXHM_Ecc()` is the eccentric-only model: the last parameter is the
initial eccentricity `e0`, defined at the reference GW frequency `f0`
(keyword argument, default 10 Hz).

This only reproduces the dominant (ℓ=2, "circular-like") harmonic's dephasing.
It does NOT include the other 9 eccentric harmonics (ℓ=1,3,4,...,10) or their
independent amplitudes (Eq. 4.29-4.31), which are a structurally separate,
much larger addition. The leading-order amplitude modulation of the dominant
harmonic (`_ecc_amp_mod` below) is defined but currently not applied in
`hphc` -- undecided whether to include it.
"""

function PolAbs(
    model::PhenomXHM_Ecc,
    f::AbstractVector,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    e0,
    Lambda1=0.0,
    Lambda2=0.0;
    f0=10.0,
)
    hp, hc = hphc(model, f, mc, eta, chi1, chi2, dL, iota, e0; f0=f0)
    return [abs.(hp), abs.(hc)]
end

function Pol(
    model::PhenomXHM_Ecc,
    f::AbstractVector,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    e0,
    Lambda1=0.0,
    Lambda2=0.0;
    f0=10.0,
)
    hp, hc = hphc(model, f, mc, eta, chi1, chi2, dL, iota, e0; f0=f0)
    return [hp, hc]
end

function Phi(
    model::PhenomXHM_Ecc,
    f,
    mc,
    eta,
    chi1,
    chi2,
    e0,
    optional_param...;
    fInsJoin_PHI=0.018,
    fcutPar=0.2,
    GMsun_over_c3=uc.GMsun_over_c3,
)
    return f .* 0.0
end

function Ampl(
    model::PhenomXHM_Ecc,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    e0;
    fcutPar=0.2,
    fInsJoin_Ampl=0.014,
    GMsun_over_c3=uc.GMsun_over_c3,
    GMsun_over_c2_Gpc=uc.GMsun_over_c2_Gpc,
    f0=10.0,
)
    hp, hc = hphc(model, f, mc, eta, chi1, chi2, dL, 0.0, e0; f0=f0)
    return abs.(hp)
end

function _fcut(
    model::PhenomXHM_Ecc,
    mc,
    eta,
    optional_param...;
    GMsun_over_c3=uc.GMsun_over_c3,
)
    return _fcut(PhenomXHM(), mc, eta; GMsun_over_c3=GMsun_over_c3)
end

function hphc(
    model::PhenomXHM_Ecc,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    e0;
    f0=10.0,
    kwargs...,
)
    hp, hc = hphc(
        PhenomXHM(),
        f,
        mc,
        eta,
        chi1,
        chi2,
        dL,
        iota;
        e0=e0,
        f0=f0,
        kwargs...,
    )

    # amp_mod = _ecc_amp_mod(e0, f, f0)  # leading-order amplitude modulation; not applied for now, see module docstring

    return hp, hc
end

function _ecc_amp_mod(e0, f, f0)

    "Leading-order-in-e0^2 amplitude modulation of the dominant harmonic."

    return 1.0 .- (157.0 / 48.0) .* e0^2 .* (f ./ f0) .^ (-19.0 / 9.0)
end
