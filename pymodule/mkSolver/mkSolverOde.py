# -*- coding: utf8 -*-
import scipy.integrate as integrate
from . import mkBaseClass
import numpy as np
"""
ODE Solver for microkinetics.
Give initial coverage, and get its time evolution till steady state
by Backward Differentiation Formula(BDF) method utilizing Scipy:
https://scipy.github.io/devdocs/reference/generated/scipy.integrate.BDF.html
"""


class MicroKineticsOde():
    def __init__(self, species, reaction, update=None):

        if update is not None and callable(update):
            self.Update = update
        elif update is not None and callable(update) is False:
            print("update para must be a function or method")

        self.species = species
        self.reaction = reaction
        self.network = {spe.name: [] for spe in self.species.values()}
        for reac in reaction:
            for spe in set(reac[0]):
                self.network[spe.name].append((-reac[0].count(spe), reac))
            for spe in set(reac[1]):
                self.network[spe.name].append((reac[1].count(spe), reac))

        self.gasSpecies = [s for s in species.values() if s.adsSpe == 0]
        self.adsSpecies = [s for s in species.values() if s.adsSpe == 1]
        # self.MakeSpeIndex()

    def Main(
        self,
        odeTime=1e3,
        odeAbsTol=1e-30,
        odeRelTol=1e-3,
    ):

        pastAdsC = [float(s.c) for s in self.adsSpecies]
        y0, t0 = [s.c for s in self.adsSpecies], 0
        r = integrate.solve_ivp(
            self.Func,
            (t0, odeTime),
            y0,
            method="BDF",
            # jac=self.Jacb,
            atol=odeAbsTol,
            rtol=odeRelTol)

        index = -1
        while min(r.y[:, index]) < 0:
            index -= 1
            if index == -len(r.t):
                msg = "Numerical explosion. Negative coverage exist in all steps of ODE."
                raise mkBaseClass.MkSolverError(msg)

        self.rall = r
        r = r.y[:, index]

        for s in self.adsSpecies:
            s.c = r[s.index]

        curAdsC = [s.c for s in self.adsSpecies]
        logErr = max(
            [abs(np.log(c1 / c2)) for c1, c2 in zip(pastAdsC, curAdsC)])
        if mkBaseClass.output >= 2:
            print(self.rall.message)
            print(self.adsSpecies)
            print(logErr)
        if mkBaseClass.output >= 1 and (index < -100):
            print("ODE Warning! negative coverage for at least 100 steps" +
                  str(index))

        if mkBaseClass.output >= 1:
            # print(rall)
            print(self.rall.message)
            print(self.Func(0, r))

        pass

    def Func(self, t, sita):
        self.Update(sita, self)
        # d\sita_i / dt = 0
        # \sigma{\sita} = 1
        func = []
        for spe in self.adsSpecies:
            f = 0
            for coeff, reac in self.network[spe.name]:
                f += reac.GenAbsRate(sita) * coeff
            func.append(f)
        # print(sita)
        # print(func)
        return func

    def PrintAll(self):
        print("index\tname\tcoverage")
        for s in self.adsSpecies.values():
            print("%5s\t%8s  %.6e" % (s.index, s.name, s.c))

    def Update(self, sita, solver):
        pass
