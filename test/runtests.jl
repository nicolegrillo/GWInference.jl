# using Pkg
# Pkg.activate("../.")
using GWInference
using Test

@testset "GWInference.jl" begin
    mc = 10.
    eta = 0.20
    chi1 = 0.1

    chi2 = 0.2
    dL = 10.
    iota = 1.

    theta = .5
    phi = 1.
    psi = 0.5

    tcoal = 0.3
    phiCoal = 0.5

    # reference values below computed with the Jan 2026 10 km ET curves (ETS), on Julia 1.10
    CE_1= CE1Id
    CE_2 = CE2NM
    ET = ETS
    network = [CE_1, CE_2, ET];
    snrD_network = SNR(PhenomD(), network, mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal) 
    @test isapprox(snrD_network, 39.156938855608196, rtol = 1e-12 )             

    snrHM_network = SNR(PhenomHM(), network, mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal)
    @test isapprox(snrHM_network, 39.3804482334247, rtol = 1e-12 )

    fisherD_network = FisherMatrix(
    PhenomD(), network, mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, phiCoal, coordinate_shift=false)
    cov = CovMatrix(fisherD_network)
    errors = Errors(cov)

    errors_tabulated= [ 0.003863803861676538
    0.0352320058894503
    0.6380118491347128
    1.6488131983913583
    0.7599165204162439
    0.031222495948098748
    0.009927687954261081
    0.06814555417938505
    0.08105308080024505
    0.008026468912306394
    1.758869875790778]

    @test isapprox(errors, errors_tabulated, rtol = 1e-8)


    fisherHM_network = FisherMatrix(
        PhenomHM(), network, mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, phiCoal, coordinate_shift=false)
    covHM = CovMatrix(fisherHM_network)
    errorsHM = Errors(covHM)

    errors_tabulatedHM= [0.001358919968470554
    0.014311145689930372
    0.27704103196699753
    0.7363589820753051
    0.7431691635712657
    0.030007509839817426
    0.009287802246134305
    0.06664432087699251
    0.08177991690258951
    0.0036155116993000025
    0.40857162217027904]

    @test isapprox(errorsHM, errors_tabulatedHM, rtol = 1e-8)
end

include("test_phenomxphm.jl")


@testset "PhenomXE" begin
    # Eccentric aligned-spin 22-mode model (port of phenomxpy IMRPhenomXE). The dynamics are integrated with an
    # adaptive ODE solver whose roundoff-level error estimate makes the accepted-step sequence (hence the waveform
    # at the ~1e-4 level) machine dependent, so the tolerances below are looser than for the other models.
    mc = 10.
    eta = 0.20
    chi1 = 0.1
    chi2 = 0.2
    dL = 10.
    iota = 1.
    theta = .5
    phi = 1.
    psi = 0.5
    tcoal = 0.3
    phiCoal = 0.5
    ecc = 0.1
    meanAno = 1.0
    network = [CE1Id, CE2NM, ETS]

    @test _npar(PhenomXE()) == 13
    @test _available_waveforms("PhenomXE") isa PhenomXE

    f = [10.0, 20.0, 50.0]
    hp, hc = hphc(PhenomXE(), f, mc, eta, chi1, chi2, dL, iota, ecc, meanAno)
    hp_tabulated = [-1.5089773722822782e-24 - 1.5704287766696464e-24im, -2.964057299671978e-25 - 8.426930779708367e-25im, 2.0315512002982394e-25 + 2.7670630692530225e-25im]
    hc_tabulated = [-1.313551870884822e-24 + 1.2621521459175385e-24im, -7.048527673427448e-25 + 2.47922291620931e-25im, 2.3144512667308263e-25 - 1.6992479503649662e-25im]
    @test isapprox(hp, hp_tabulated, rtol = 1e-3)
    @test isapprox(hc, hc_tabulated, rtol = 1e-3)
    # the e -> 0 waveform has the same amplitude as IMRPhenomXAS
    hp0, hc0 = hphc(PhenomXE(), f, mc, eta, chi1, chi2, dL, iota, 0.0, 0.0)
    ampXAS = Ampl(PhenomXAS(), f, mc, eta, chi1, chi2, dL)
    @test isapprox(abs.(hp0), 0.5 .* (1 .+ cos(iota)^2) .* ampXAS, rtol = 1e-3)

    snrXE_network = SNR(PhenomXE(), network, mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, ecc, meanAno)
    @test isapprox(snrXE_network, 39.0370542106737, rtol = 1e-3)

    fisherXE_network = FisherMatrix(
        PhenomXE(), network, mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, phiCoal, ecc, meanAno, coordinate_shift=false)
    @test size(fisherXE_network) == (13, 13)
    @test all(isfinite, fisherXE_network)
    errorsXE = Errors(CovMatrix(fisherXE_network))

    errors_tabulatedXE = [0.0012796717246793979
    0.0018435709895120481
    0.08938538597397967
    0.31084094216970887
    0.8498699573206085
    0.03204263991447796
    0.009452124375376394
    0.07617285762382942
    0.0799775690624142
    0.00016269317307005806
    7.445299457352858
    0.00045052378284311096
    2.7952757253895717]

    @test isapprox(errorsXE, errors_tabulatedXE, rtol = 1e-1)
end