# -*- coding: utf8 -*-
import numpy as np
import mpmath as mp
# temperature
temp = 500
# flag to lable adsorbates. Support only one site currently
adsflag = ["ads", "*", "$"]
# 0 to ignore most output
output = -1
#
energyType = "energy"


def KbT():
    """k_bT in eV"""
    return 8.314462618 * temp / 96485


def KbTh():
    """k_bT / h"""
    return KbT() / 4.1356676969e-15


class Spec():
    def __init__(self, name, **kargs):
        self.name = name
        if "energy" in kargs:
            self.trueEnergy = kargs["energy"]  # energy / enthapy
            self.energy = kargs["energy"]  # gibbs energy in fact
        else:
            # pass
            self.energy = 0
            self.trueEnergy = 0

        if "entropy" in kargs:
            self.entropy = kargs["entropy"]
        else:
            self.entropy = 0
        self.c = 0  # coverage
        self.index = None

    def __repr__(self):
        return "Spec_%s_%.5e" % (self.name, self.c)

    def CalT(self):
        if energyType == "energy":
            self.energy = self.trueEnergy - self.entropy * temp
        else:
            pass
            # self.energy = self.trueEnergy - self.entropy * temp


class Reac(list):
    def __init__(self, species1, species2, **kargs):
        list.__init__(self)
        self.append(species1)
        self.append(species2)
        if "tsEnergy" in kargs:
            self.trueTS = kargs["tsEnergy"]
            self.TS = kargs["tsEnergy"]
        elif "bf" in kargs and "bb" in kargs:
            self.bf = kargs["bf"]
            self.bb = kargs["bb"]
        if "entropy" in kargs:
            self.entropy = kargs["entropy"]
        else:
            self.entropy = 0

        self.fCount = {s: self[0].count(s) for s in set(self[0])}
        self.bCount = {s: self[1].count(s) for s in set(self[1])}
        pass

    def __repr__(self):
        return "Rect_" + ".".join(s.name for s in self[0]) + "_" + ".".join(
            s.name for s in self[1])

    def CalT(self):
        if energyType == "energy":
            # TS is forward Barrier
            self.TS = self.trueTS - self.entropy * temp
            # forward barrier
            self.bf = self.TS - sum([s.energy for s in self[0]])
            # backward barrier
            self.bb = self.TS - sum([s.energy for s in self[1]])
            self.kf = np.exp(-self.bf / KbT()) * KbTh()
            self.kb = np.exp(-self.bb / KbT()) * KbTh()
        elif energyType == "barrier":
            dSfT = temp * (self.entropy - sum([s.entropy for s in self[0]]))
            dsbT = temp * (self.entropy - sum([s.entropy for s in self[1]]))
            self.gf = self.bf - dSfT
            self.gb = self.bb - dsbT
            self.kf = np.exp(-(self.bf - dSfT) / KbT()) * KbTh()
            self.kb = np.exp(-(self.bb - dsbT) / KbT()) * KbTh()

    # def GenRate(self, flag=0):
    #     """
    #     return r_+ for flag=0
    #     retrun r_- for flag=1
    #     """
    #     if flag == 0:
    #         return KbTh() * self.kf * np.prod([s.c for s in self[0]])
    #     elif flag == 1:
    #         return KbTh() * self.kb * np.prod([s.c for s in self[1]])
    #     else:
    #         raise BaseException

    def GenAbsRate(self, sita):
        rabs = self.kf * np.prod([sita[s.index] if s.index is not None else s.c for s in self[0]]) - \
             self.kb * np.prod([sita[s.index] if s.index is not None else s.c for s in self[1]])
        return rabs

    def GenRatePartial(self, sita, spe):
        if spe in self[0]:
            n = self.fCount[spe]
            drds = self.kf * np.prod([sita[s.index] if s.index is not None else s.c for s in self[0] if s is not spe]) * \
                n * sita[spe.index] ** (n-1)
        elif spe in self[1]:
            n = self.bCount[spe]
            drds = self.kb * np.prod([sita[s.index] if s.index is not None else s.c for s in self[1] if s is not spe]) * \
                n * sita[spe.index] ** (n-1)
        else:
            drds = 0
        return drds

    def GenBarrierHeight(self):
        muf = self.bf - KbT() * np.sum([float(mp.log(s.c)) for s in self[0]])
        mub = self.bb - KbT() * np.sum([float(mp.log(s.c)) for s in self[1]])
        return (muf, mub)


class MkSolverError(Exception):
    def __init__(self, msg):
        self.message = msg

    def __repr__(self):
        # print(self.message)
        return self.message

    def __str__(self):
        return self.__repr__()
