# -*- coding: utf8 -*-
from . import mkBaseClass
from .mkSinglePoint import ReadReactInput, Output, SingleMK
from . import mkOde, mkOpt


class SingleMKOrga(SingleMK):
    def __init__(
        self,
        species,
        reaction,
        temperature,
        concentration,
    ):
        super().__init__(species,
                         reaction,
                         temperature,
                         concentration,
                         update=False,
                         initC=True)

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

    def Run(
        self,
        odeTime=1e4,
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
            print(str(self.temp) + " " + str(self.gasPressure) + " start")
        # mk1 = mkOde(self.species, self.reaction, update=self.Update)
        mk1 = mkOde(self.species, self.reaction)
        if odeTime > 0:
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

        mk2 = mkOpt(species, reaction)
        error = float(min(mk2.Func([s.c for s in mk2.adsSpecies])))
        if mkBaseClass.output >= 0:
            print(
                str(self.temp) + " " + str(self.gasPressure) +
                "Error to Steady state: Err=%5e" % error)
        self.error = float(error)
        adsRate = {
            spe.name: sum([
                reac.GenAbsRate([s.c for s in mk2.adsSpecies]) * coeff
                for coeff, reac in mk2.network[spe.name]
            ])
            for spe in mk2.adsSpecies
        }
        adsCvrg = {spe.name: spe.c for spe in mk2.adsSpecies}
        self.result = (self.temp, self.gasPressure, adsRate, adsCvrg)
        minC = min(adsCvrg.values())
        if minC > 0:
            self.success = True
        else:
            self.success = False
        return self
