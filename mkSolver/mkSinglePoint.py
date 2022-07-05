# -*- coding: utf8 -*-
from . import Spec, Reac, mkOde, mkOpt, mkBaseClass
import numpy as np
import scipy.optimize as opt


def ReadReactInput(filename, energyType="energy"):
    """
    Read reaction file to get relative energy of the microkinetics model
    support two input type:
    1. energyType = "energy", that input energy are the height on energy profile of each species.
    2. energyType = "barrier", that input are forward and backward energy barrier
    """
    mkBaseClass.energyType = energyType
    flag = "species"
    species = {}
    reaction = []

    with open(filename) as fp:
        lines = fp.readlines()

    for line in lines:
        if line.startswith("#"):
            continue
        li = line.lower().split()
        if len(li) == 0:
            continue
        if li[0] == "species":
            flag = li[0]
            continue
        if li[0] == "reaction":
            flag = li[0]
            continue

        if flag == "species":
            if energyType == "energy":
                if len(li) == 2:
                    species[li[0]] = Spec(li[0], energy=float(li[1]))
                elif len(li) == 3:
                    species[li[0]] = Spec(li[0],
                                          energy=float(li[1]),
                                          entropy=float(li[2]))
                else:
                    print("Wrong species Input")
            elif energyType == "barrier":
                if len(li) == 1:
                    species[li[0]] = Spec(li[0])
                elif len(li) == 2:
                    species[li[0]] = Spec(li[0], entropy=float(li[1]))
                else:
                    print("Wrong species Input")
        elif flag == "reaction":
            if energyType == "energy":
                dash = li.index("-")
                try:
                    float(li[-2])
                    numb = -2
                except BaseException:
                    numb = -1

                compos1 = li[0:dash]
                spe1 = [species[s] for s in compos1]
                compos2 = li[dash + 1:numb]
                spe2 = [species[s] for s in compos2]
                if numb == -1:
                    pair = Reac(spe1, spe2, tsEnergy=float(li[-1]))
                elif numb == -2:
                    pair = Reac(spe1,
                                spe2,
                                tsEnergy=float(li[-2]),
                                entropy=float(li[-1]))
                else:
                    print("Wrong reaction Input")
                reaction.append(pair)
            elif energyType == "barrier":
                dash = li.index("-")
                try:
                    float(li[-3])
                    numb = -3
                except BaseException:
                    numb = -2

                compos1 = li[0:dash]
                spe1 = [species[s] for s in compos1]
                compos2 = li[dash + 1:numb]
                spe2 = [species[s] for s in compos2]
                if numb == -2:
                    pair = Reac(spe1, spe2, bf=float(li[-2]), bb=float(li[-1]))
                elif numb == -3:
                    pair = Reac(spe1,
                                spe2,
                                bf=float(li[-3]),
                                bb=float(li[-2]),
                                entropy=float(li[-1]))
                else:
                    print("Wrong reaction Input")
                reaction.append(pair)

    return (species, reaction)


def Output(tasks, outFile="out.txt"):
    col = ["temp"]
    col.extend([
        spe.name + "Prsr" for spe in tasks[0].species.values()
        if spe.adsSpe == 0
    ])
    col.extend([
        spe.name + "Rate" for spe in tasks[0].species.values()
        if spe.adsSpe == 0
    ])
    col.extend([
        spe.name + "Cvrg" for spe in tasks[0].species.values()
        if spe.adsSpe == 1
    ])
    col = {c: i for i, c in enumerate(col)}
    table = [[] for k in col]
    result = [t.result for t in tasks]
    for r in result:
        if r is None:
            continue
        table[col["temp"]].append(r[0])
        for k, v in r[1].items():
            table[col[k + "Prsr"]].append(v)
        for k, v in r[2].items():
            table[col[k + "Rate"]].append(v)
        for k, v in r[3].items():
            table[col[k + "Cvrg"]].append(v)

    out = ["\t".join(["%12s" % c for c in col]) + "\n"]
    for i in range(len(table[0])):
        out.append("\t".join(["%8e" % table[col[c]][i] for c in col]) + "\n")
    print("\n".join(out))
    with open(outFile, "w") as fp:
        fp.writelines(out)


class SingleMK():
    def __init__(self,
                 species,
                 reaction,
                 temperature,
                 gasPressure,
                 update=None,
                 initC=True):
        """
        species: (dict) { (str) speciesName: (SpecObject) species  }
        reaction (list) [ (ReacObject) reaction]
        temperature: (float) temp
        gasPressure: (dict) { (str) gasSpeciesName: (float) pressureInAtm}
        """
        self.species = species
        self.reaction = reaction
        self.gasPressure = gasPressure
        self.temp = temperature
        self.containSite = {}  # existing adsflag, e.g., "#", "ads"
        self.initC = initC

        if initC:
            self.InitSpecies()
        self.gasSpecies = [s for s in species.values() if s.adsSpe == 0]
        self.adsSpecies = [s for s in species.values() if s.adsSpe == 1]
        # if calT:
        if update and callable(update) and not initC:
            self.Update = update
        else:
            self.Update = None

    def Run(
        self,
        odeTime=1e7,
        opt=True,
        odeTail=True,
        odeAbsTol=1e-30,
        odeRelTol=1e-3,
    ):
        if self.initC:
            if self.temp is not None:
                mkBaseClass.temp = self.temp
            if self.gasPressure is not None:
                for speName in self.gasPressure:
                    self.species[speName].c = self.gasPressure[speName]
            for spe in self.species.values():
                spe.CalT()
            for reac in self.reaction:
                reac.CalT()
        if self.Update and callable(self.Update) and not self.initC:
            self.Update(self.species, self.reaction)

        if mkBaseClass.output >= 0:
            print(
                str(self.temp) + " " + str(self.gasPressure) + " start")
        # mk1 = mkOde(self.species, self.reaction, update=self.Update)
        mk1 = mkOde(self.species, self.reaction)
        if odeTime > 0:
            self.CoverageNormal()
            mk1.Main(odeTime=odeTime, odeAbsTol=odeAbsTol, odeRelTol=odeRelTol)
            species, reaction = mk1.species, mk1.reaction
            self.odeResult = mk1.rall
            if mkBaseClass.output >= 1:
                print("Ode Ouput")
                print([(spe.name,
                        sum([
                            reac.GenAbsRate([s.c
                                             for s in mk1.adsSpecies]) * coeff
                            for coeff, reac in mk1.network[spe.name]
                        ])) for spe in mk1.gasSpecies])
                print(mk1.adsSpecies)

            if mkBaseClass.output >= 0:
                print(
                    str(self.temp) + " " + str(self.gasPressure) + " ODE Done")

        if odeTail:
            self.CoverageNormal()
            mk1.Main(odeTime=10, odeAbsTol=odeAbsTol, odeRelTol=odeRelTol)
            species, reaction = mk1.species, mk1.reaction
            if mkBaseClass.output >= 1:
                print("Ode Ouput")
                print([(spe.name,
                        sum([
                            reac.GenAbsRate([s.c
                                             for s in mk1.adsSpecies]) * coeff
                            for coeff, reac in mk1.network[spe.name]
                        ])) for spe in mk1.gasSpecies])
                print(mk1.adsSpecies)
            if mkBaseClass.output >= 0:
                print(
                    str(self.temp) + " " + str(self.gasPressure) + " ODE Done")

        mk2 = mkOpt(species, reaction)
        # mk2 = mkOpt(species, reaction)
        if opt:
            self.CoverageNormal()
            mk2.Main()
            error = float(min(mk2.Func([s.c for s in mk2.adsSpecies])))
            if mkBaseClass.output >= 1:
                print("Opt Ouput")
                print([(spe.name,
                        sum([
                            reac.GenAbsRate([s.c
                                             for s in mk2.adsSpecies]) * coeff
                            for coeff, reac in mk2.network[spe.name]
                        ])) for spe in mk2.gasSpecies])
            if mkBaseClass.output >= 0:
                print(
                    str(self.temp) + " " + str(self.gasPressure) +
                    " Optimization Done, Err=%5e" % error)
        else:
            error = float(min(mk2.Func([s.c for s in mk2.adsSpecies])))
            if mkBaseClass.output >= 0:
                print(
                    str(self.temp) + " " + str(self.gasPressure) +
                    " Skip Optimization, Err=%5e" % error)
        self.error = float(error)
        gasRate = {
            spe.name: sum([
                reac.GenAbsRate([s.c for s in mk2.adsSpecies]) * coeff
                for coeff, reac in mk2.network[spe.name]
            ])
            for spe in mk2.gasSpecies
        }
        adsCvrg = {spe.name: spe.c for spe in mk2.adsSpecies}
        self.result = (self.temp, self.gasPressure, gasRate, adsCvrg)
        minC = min(adsCvrg.values())
        if minC > 0:
            self.success = True
        else:
            self.success = False
        return self

    def InitSpecies(self):
        # init coverage
        # set adsSpe=1 for adsorbates, set adsSpe=0 for gas
        # set coverage=0.001 for adsorbates, set coverage=1 for gas & site
        for name, spe in self.species.items():
            flag = False
            for f in mkBaseClass.adsflag:
                if f in name:
                    spe.adsSpe = 1
                    spe.c = 0
                    flag = True
            if flag is False:
                spe.adsSpe = 0
                if hasattr(spe, "c") is False or spe.c == 0:
                    spe.c = 1
        for f in mkBaseClass.adsflag:
            if f in self.species:
                self.species[f].c = 1
                self.species[f].site = True
                if not hasattr(self.species[f], "total"):
                    self.species[f].total = 1

        # add species index
        index = 0
        for spe in self.species.values():
            if spe.adsSpe == 1:
                # adsorbates
                spe.index = index
                index += 1
            else:
                # gas species
                spe.index = None

        for f in mkBaseClass.adsflag:
            if sum([
                    name.count(f)
                    for name, spe in self.species.items() if spe.adsSpe == 1
            ]) != 0:
                self.containSite[f] = self.species[f].total

    def CoverageNormal(self):
        # Nomallize the sum of each coverage site to 1
        allC = np.array(
            [spe.c for name, spe in self.species.items() if spe.adsSpe == 1])

        bounds = opt.Bounds([0] * len(allC), [1] * len(allC))

        r = opt.minimize(lambda X: sum([(np.dot(
            np.array([
                name.count(f) for name, spe in self.species.items()
                if spe.adsSpe == 1
            ]), X) - self.containSite[f])**2 for f in self.containSite]),
                         allC,
                         method="SLSQP",
                         bounds=bounds,
                         options={'ftol': 1e-10})
        for s in self.species.values():
            if s.adsSpe == 1:
                s.c = r.x[s.index]
        return

    def Update(self, sita, solver):
        pass
