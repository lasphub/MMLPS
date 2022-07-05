from numpy import log

pathSampleMultiPara = {
    "direcToBallance": {
        "/home10/shiyf/3.CO2Hydro/0320_test/1.CO2+H2": "H2.H2",
        "/home10/shiyf/3.CO2Hydro/0320_test/3.HCHO+H2": "H2O",
        "/home10/shiyf/3.CO2Hydro/0320_test/2.HCOOH+H2": "H2",
    },
    "direcJoint": "/home10/shiyf/3.CO2Hydro/0320_test/joint"
}

pathSamplePara = {
    "prog": "/home10/shiyf/prog/FULL_pro2.2.3_intel18_sss/Src/lasp",
    "runmode": "local",
    "runcycle": 30,
    "SSWstep": 2000,
    "ncputotal": 96,  # for local mode
    "ncpusample": 24,
    "ncpusplit": 8,
    "MaxNumofSinglePair": 50,
    "ThresholdofSelect": 500,
    "Lscreenupper": 1,
    "restart": 1,
}


GibbsCorrDict = {
    # NIST exp data
    # https://janaf.nist.gov/
    # - T * S + [ H(T) - H(0K) ] + kb * T * ln(P/P0)
    "H2":
    -500 * 145.737 * 1.036e-5 + (5.882 + 8.467) * 1.036e-2 +
    8.61734e-5 * 500 * log(40 / 1),  # + 0.27,
    "H2O":
    -500 * 206.534 * 1.036e-5 + (6.925 + 9.904) * 1.036e-2 +
    8.61734e-5 * 500 * log(1 / 1),  # + 0.52,
    "CO2":
    -500 * 234.901 * 1.036e-5 + (8.305 + 9.364) * 1.036e-2 +
    8.61734e-5 * 500 * log(10 / 1),  # + 0.31,
    "CO":
    -500 * 212.831 * 1.036e-5 + (5.931 + 8.671) * 1.036e-2 +
    8.61734e-5 * 500 * log(10 / 1) - 0.44,  # + 0.13,
    "CH3OH":
    0.202312 - 1.355580 + 8.61734e-5 * 500 * log(1 / 1),
}
