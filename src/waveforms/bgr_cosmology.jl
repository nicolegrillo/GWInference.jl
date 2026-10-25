"""
Cosmology helpers shared by the beyond-GR XHM models (EdGB, dCS, MG).

Included once from waveforms_master.jl, before PhenomXHM_EdGB.jl / PhenomXHM_dCS.jl / PhenomXHM_MG.jl.

* EdGB and dCS use `get_z(dL_Gpc)` to turn the detector-frame mass into the source-frame mass (the couplings alpha
  have units of length^2 and enter through the physical masses).
* MG uses `get_z` and `get_dL_MG`, the distance D of Will, Eq. (2.5) of gr-qc/9709011, which is NOT the luminosity distance.

The constants carry a `_bgr_` prefix so they cannot clash with the ones in utils.jl (where H0 carries Unitful units).
"""

using QuadGK

const _bgr_clight_kms = 2.99792458e5 # km/s

const _bgr_Omega0_m = 0.3153

const _bgr_Omega0_Lambda = 1 - _bgr_Omega0_m

const _bgr_H0 = 67.66 # km/s 1/Mpc

"""
Hubble constant as a function of redshift, results in km/s 1/Mpc

"""
function get_H_z(z, H0, Omega0_m, Omega0_Lambda)
    H_z = H0 * (Omega0_m * (1 .+ z) .^ 3 .+ Omega0_Lambda) .^ 0.5 # km/s 1/Mpc
    return H_z
end

"""
Luminosity distance as a function of redshift, results in Mpc. clight in km/s

"""
function get_dL(z, clight = _bgr_clight_kms, H0 = _bgr_H0, Omega0_m = _bgr_Omega0_m, Omega0_Lambda = _bgr_Omega0_Lambda)
    dL = zeros(length(z))
    jj = 1
    for ii in z
        dL[jj] =
            quadgk(k -> 1.0 ./ get_H_z(k, H0, Omega0_m, Omega0_Lambda), 0.0, ii)[1] *
            clight *
            (1.0 + ii)
        jj += 1
    end
    return dL
end

"""
Distance D entering the massive-graviton phase (Will 1997, Eq. 2.5): D = (1+z) int dz' / ((1+z')^2 H(z')), results in Mpc.
It is NOT the luminosity distance. clight in km/s

"""
function get_dL_MG(z, clight = _bgr_clight_kms, H0 = _bgr_H0, Omega0_m = _bgr_Omega0_m, Omega0_Lambda = _bgr_Omega0_Lambda)
    dL = Vector{promote_type(eltype(z), Float64)}(undef, length(z))
    jj = 1
    for ii in z
        z0 = ForwardDiff.value(ii)   # quadgk cannot integrate up to a Dual: evaluate at the value, add the analytic derivative
        I = quadgk(k -> 1.0 / ((1.0 + k)^2 * get_H_z(k, H0, Omega0_m, Omega0_Lambda)), 0.0, z0)[1]
        val = I * clight * (1.0 + z0)
        slope = clight * (I + 1.0 / ((1.0 + z0) * get_H_z(z0, H0, Omega0_m, Omega0_Lambda)))   # d/dz of (1+z) int dz'/((1+z')^2 H)
        dL[jj] = val + slope * (ii - z0)
        jj += 1
    end
    return dL
end

"""
Redshift as a function of luminosity distance (in Gpc). clight in km/s
"""
function get_z(dL)

    # use bisection method to find z
    z=1
    zmin = 0.
    zmax = 30.
    dL = dL*1e3 # convert to Mpc
    iteration = 0

        while abs.(dL - get_dL(z)[1]) .> 1e-3
            iteration += 1

            if iteration > 100
                println("Could not find z")
                break
            end

            z = 0.5*(zmin + zmax)
            if dL .> get_dL(z)[1]
                zmin = z
            else
                zmax = z
            end
        end
    # one Newton step so that z carries the derivative wrt dL (value unchanged to working precision)
    slope = _bgr_clight_kms * (quadgk(k -> 1.0 / get_H_z(k, _bgr_H0, _bgr_Omega0_m, _bgr_Omega0_Lambda), 0.0, z)[1] + (1.0 + z) / get_H_z(z, _bgr_H0, _bgr_Omega0_m, _bgr_Omega0_Lambda))
    z = z + (dL - get_dL(z)[1]) / slope
    return z
end
