# -*- coding: utf8 -*-
import scipy.optimize as opt
from mpmath import mp
import mpmath
from . import mkBaseClass

mp.dps = 50
"""
Optimization Solver for microkinetics.
Assume at steady state, and optimize coverage to minimize function errors
from well initial guess by MINPACK’s hybrd and hybrj routines (modified Powell method) method, utilizing Scipy,
combining with mpmath MDNewton optimization for high precision solution
"""


class MicroKineticsOpt():
    def __init__(self, species, reaction, update=None):

        self.cycleMax = 1e6
        self.cycle = 1
        self.minCoverage = 1e-100

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

    def Main(self):
        r = opt.root(self.Func, [s.c for s in self.adsSpecies if s.name],
                     method='hybr',
                     options={"xtol": 1e-12})

        for s in self.adsSpecies:
            s.c = r.x[s.index]
        if mkBaseClass.output >= 2:
            print("ScipyOpt: " + r.message)
            print(r.fun)

        try:
            r = mpmath.findroot(self.Func2,
                                [s.c for s in self.adsSpecies if s.name],
                                verify=False)
            r = list(r[:, 0])

            error = float(min([float(x) for x in self.Func(r)]))
            if abs(error) > 1e-10:
                gasPressure = {spe.name: spe.c for spe in self.gasSpecies}
                print(
                    str(mkBaseClass.temp) + " " + str(gasPressure) +
                    " Warning: Error of %.8e > 1e-10" % error)

            for s in self.adsSpecies:
                s.c = r[s.index]

            if mkBaseClass.output >= 2:
                print("mpmathOpt:")
                # print([float(x) for x in self.Func(r)])
                print(self.adsSpecies)

        except ZeroDivisionError as e:
            print("mpmath Error:")
            print(e)
        pass

    def Func(self, sita):
        self.Update(sita, self)
        func = []
        for spe in self.adsSpecies:
            f = 0
            if hasattr(spe, "site") and spe.site is True:
                # if spe is a site, calculate the sum of site
                for _s in self.adsSpecies:
                    f += _s.name.count(spe.name) * sita[_s.index]
                f -= spe.total
            else:
                # if spe is an adsotbates, calculate it total reaction rate
                for coeff, reac in self.network[spe.name]:
                    f += reac.GenAbsRate(sita) * coeff
            func.append(f)
        return func

    def Func2(self, *sita):
        self.Update(sita, self)
        func = []
        for spe in self.adsSpecies:
            f = 0
            if hasattr(spe, "site") and spe.site is True:
                # if spe is a site, calculate the sum of site
                for _s in self.adsSpecies:
                    f += _s.name.count(spe.name) * sita[_s.index]
                f -= spe.total
            else:
                # if spe is an adsotbates, calculate it total reaction rate
                for coeff, reac in self.network[spe.name]:
                    f += reac.GenAbsRate(sita) * coeff
            func.append(f)
        return func

    def Update(self, sita, solver):
        pass
