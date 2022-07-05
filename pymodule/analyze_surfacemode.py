#! /usr/bin/env

import os
# from allstr_new import allstr
import numpy as np
import re
import hashlib
# import matplotlib.pyplot as plt
# from reactpair import Pair as allpair
# import pickle
# from ECFP import ECFP
from ECFP import AllStr as allstr_ecfp
import ctypes
from ctypes import pointer
from math import log
import networkx as nx
import copy as cp
from multiprocessing import Pool
from pathsample_para import GibbsCorrDict
import sys

class Pair(list):
    def __init__(self, str1, str2, TS, TSstr, sortflag=True):
        list.__init__(self)
        self.append(str1)
        self.append(str2)
        # self.append(TSstr)
        self.TS = TS
        if sortflag:
            self.sort(key=lambda X: X.name)
        self.name = '%s' % (self[0].name + self[1].name)
        self.TSstr = TSstr
        return

    def __hash__(self):
        return id(self)

    def ISFSdistance(
            self,
            program='/home10/kpl/pymodule/reactioncoord/reactioncoord.so'):
        self[0].calCtypes()
        self[1].calCtypes()

        self.factxyz = pointer(
            (ctypes.c_double * 9)(*[0, 0, 0, 0, 0, 0, 0, 0, 0]))
        self.length = pointer((ctypes.c_double)(0))
        len3n3 = self[0].natom * 3 + 9
        self.fracIS = pointer((ctypes.c_double * len3n3)(*range(len3n3)))
        self.fracFS = pointer((ctypes.c_double * len3n3)(*range(len3n3)))

        #   print self[0].c_rv.contents
        #   for item in self[0].c_rv.contents:
        #       print item
        #   print self[0].c_xa.contents
        #   for item in self[0].c_xa.contents:
        #       print item
        self.rcal = ctypes.cdll.LoadLibrary(program)
        self.rcal.isfslength_new_(self[0].c_natm, self[0].c_rv, self[0].c_xa,
                                  self[1].c_rv, self[1].c_xa, self.factxyz,
                                  self.length, self.fracIS, self.fracFS)
        # print self.length.contents
        # print self.length[0]
        self.totdis = self.length[0]
        return self.length[0]

    def ISFSdistance_crystal(
            self,
            program='/home10/kpl/pymodule/reactioncoord/reactioncoord.so'):
        self[0].calCtypes()
        self[1].calCtypes()

        self.factxyz = pointer(
            (ctypes.c_double * 9)(*[0, 0, 0, 0, 0, 0, 0, 0, 0]))
        self.length = pointer((ctypes.c_double)(0))
        len3n3 = self[0].natom * 3 + 9
        self.fracIS = pointer((ctypes.c_double * len3n3)(*range(len3n3)))
        self.fracFS = pointer((ctypes.c_double * len3n3)(*range(len3n3)))

        self.rcal = ctypes.cdll.LoadLibrary(program)
        self.rcal.isfslength_crystal_(self[0].c_natm, self[0].c_rv,
                                      self[0].c_xa, self[1].c_rv, self[1].c_xa,
                                      self.factxyz, self.length, self.fracIS,
                                      self.fracFS)
        # print self.length.contents
        # print self.length[0]
        self.totdis = self.length[0]
        return self.length[0]

    def GetFS(self, nameIS):
        if self[0].name == nameIS:
            return self[1]
        elif self[1].name == nameIS:
            return self[0]
        else:
            # print("Err: Wrong Reactang")
            return None


class Path(list):
    def __init__(self):
        list.__init__(self)
        self.maxTS = 999

    def __hash__(self):
        dim = 10000000
        md5 = hashlib.md5()
        md5.update(self[0][0].name.encode())
        for pair in self:
            md5.update(pair[1].name.encode())
        _tmp = md5.hexdigest()
        return int((int(_tmp, 16)) % dim)

    def __repr__(self):
        name = [self[0][0].sminame]
        for pair in self:
            name.append(pair[1].sminame)
        return "\n".join(["%2d  " % i + n for i, n in enumerate(name)])

    def copy(self, cpath):
        for pair in cpath:
            self.append(pair)
        self.maxTS = cpath.maxTS
        self.namechain = cpath.namechain[:]

    def printpath(self, fname):
        _path = allstr()
        for pair in self:
            _path.append(pair[0])
            _path.append(pair[1])
        _path.printall(fname)

    def printpath_TS(self, fname):
        _path = allstr()
        _path.append(self.allstrmin[self[0][0].name])
        for pair in self:
            _path.append(pair[0])
            _path.append(pair.TSstr)
            _path.append(pair[1])
            _path.append(self.allstrmin[pair[1].name])
        _path.printall(fname)

    def GetBarrier(self):
        Ecurrent = 0
        bmax = 0
        for pair in self:
            if pair[0].energy < Ecurrent:
                Ecurrent = pair[0].energy
            bcurrent = pair.TS - Ecurrent
            if bmax < bcurrent:
                bmax = bcurrent
            if pair[1].energy < Ecurrent:
                Ecurrent = pair[1].energy
        self.barrier = bmax
        self.score = bmax * 1000 + len(self)
        return

    def GetBarrier_new(self, allmin):
        Ecurrent = 0
        bmax = 0
        barrierpair = 0
        for ipair, pair in enumerate(self):
            if allmin[pair[0].name].energy < Ecurrent:
                Ecurrent = allmin[pair[0].name].energy
            bcurrent = pair.TS - Ecurrent
            if bmax < bcurrent:
                bmax = bcurrent
                barrierpair = ipair
            if allmin[pair[1].name].energy < Ecurrent:
                Ecurrent = allmin[pair[1].name].energy
        self.barrier = bmax
        self.score = bmax * 1000 + len(self)
        self.barrierpair = barrierpair
        return

    def GetBarrier_list(self, allmin):

        self.GetBarrier_new(allmin)
        barrierpair = self.barrierpair
        self.barrierlist = [self.barrier]

        while (barrierpair + 1) != len(self):
            _tmppath = Path()
            for ipair in range(barrierpair + 1, len(self)):
                _tmppath.append(self[ipair])
            _tmppath.GetBarrier_new(allmin)
            barrierpair = barrierpair + 1 + _tmppath.barrierpair
            self.barrierlist.append(_tmppath.barrier)

        # print self.barrierlist
        self.score = 0
        if len(self.barrierlist) > 15:
            print('maybe too long to cal score')
        for i in range(min(15, len(self.barrierlist))):
            self.score = self.score + (self.barrierlist[i] // 0.01 /
                                       100) * (1000**(4 - i))
        self.score = self.score + len(self)

        return

    def EndtoendSort(self, name1):
        for i, pair in enumerate(self):
            pair.sort(key=lambda X: X.name)
            if pair[0].name == name1:
                name1 = pair[1].name
            else:
                # print 'reverse', i
                # print pair[0].name,pair[1].name
                pair.sort(key=lambda X: X.name, reverse=True)
                # if pair[0].name != name1:
                #    pair.sort(key= lambda X: X.name)
                # print pair[0].name,pair[1].name
                name1 = pair[1].name
                self.pop(i)
                self.insert(i, pair)

    # def plot(self, fname):
    #     IS, TS, FS = [], [], []
    #     line, linex = [], []
    #     for i, pair in enumerate(self):
    #         IS.append(pair[0].energy)
    #         TS.append(pair.TSstr.energy)
    #         FS.append(pair[1].energy)
    #         line.append([pair[0].energy, pair.TSstr.energy, pair[1].energy])
    #         linex.append([i + 1, i + 1.5, i + 2])
    #     plotISx = np.linspace(1, len(self), len(self))
    #     plotTSx = np.linspace(1.5, len(self) + 0.5, len(self))
    #     plotFSx = np.linspace(2, len(self) + 1, len(self))
    #     plt.figure(figsize=(8, 6), dpi=400)
    #     plt.scatter(plotISx, IS, s=400, color='black', marker='_', lw=5)
    #     plt.scatter(plotTSx, TS, s=400, color='red', marker='_', lw=5)
    #     plt.scatter(plotFSx, FS, s=400, color='black', marker='_', lw=5)
    #     for pairlinex, pairline in zip(linex, line):
    #         plt.plot(pairlinex, pairline, 'b-')

    #     # plt.xticks([])
    #     # plt.yticks()
    #     plt.xlabel('Reaction Coordinate', fontsize=24, fontweight='bold')
    #     plt.ylabel('Delta [eV per atom]', fontsize=24, fontweight='bold')
    #     plt.savefig(fname, dpi=400)

    #     return

    # def plot_new(self, fname, allmin):
    #     IS, TS, FS, ISmin, FSmin = [], [], [], [], []
    #     line, linex = [], []
    #     for i, pair in enumerate(self):
    #         IS.append(pair[0].energy)
    #         ISmin.append(allmin[pair[0].name].energy)
    #         TS.append(pair.TSstr.energy)
    #         FS.append(pair[1].energy)
    #         FSmin.append(allmin[pair[1].name].energy)
    #         line.append([pair[0].energy, pair.TSstr.energy, pair[1].energy])
    #         linex.append([i + 1, i + 1.5, i + 2])
    #     plotISx = np.linspace(1, len(self), len(self))
    #     plotTSx = np.linspace(1.5, len(self) + 0.5, len(self))
    #     plotFSx = np.linspace(2, len(self) + 1, len(self))
    #     plt.figure(figsize=(8, 6), dpi=400)
    #     plt.scatter(plotISx, IS, s=400, color='black', marker='_', lw=5)
    #     plt.scatter(plotTSx, TS, s=400, color='red', marker='_', lw=5)
    #     plt.scatter(plotFSx, FS, s=400, color='black', marker='_', lw=5)
    #     plt.scatter(plotISx, ISmin, s=400, color='green', marker='_', lw=5)
    #     plt.scatter(plotFSx, FSmin, s=400, color='green', marker='_', lw=5)
    #     for pairlinex, pairline in zip(linex, line):
    #         plt.plot(pairlinex, pairline, 'b-')

    #     # plt.xticks([])
    #     # plt.yticks()
    #     plt.xlabel('Reaction Coordinate', fontsize=24, fontweight='bold')
    #     plt.ylabel('Delta [eV per atom]', fontsize=24, fontweight='bold')
    #     plt.savefig(fname, dpi=400)

    #     return

    def plot_paper3(self, allmin):
        point = []
        zeroe = allmin[self[0][0].name].energy
        print(zeroe)
        point = [0, 0]
        for i, pair in enumerate(self):
            # if pair.TSstr.energy-zeroe > 3000 : zeroe = allmin[pair[0].name].energy-point[-1]
            point.append(pair.TSstr.energy - zeroe)
            point.append(pair.TSstr.energy - zeroe)
            point.append(allmin[pair[1].name].energy - zeroe)
            point.append(allmin[pair[1].name].energy - zeroe)
            print(allmin[pair[1].name].energy)
            print(pair[1].name)
            print(pair.TSstr.energy - allmin[pair[0].name].energy)
            print(allmin[pair[1].name].energy - allmin[pair[0].name].energy)

        plotx = np.linspace(1 - 0.125,
                            len(self) + 1 + 0.125, 4 * len(self) + 2)
        return plotx, point

    def plot_forcecontinue(self, allmin, outfile='outpath'):
        point = []
        zeroe = allmin[self[0][0].name].energy
        # print(zeroe)
        point = [0, 0]
        FS = 0

        tmp = allstr()
        # tmp.append(allmin[self[0][0].name])

        for i, pair in enumerate(self):
            zeroe = allmin[pair[0].name].energy - FS
            # if pair.TSstr.energy-zeroe > 3000 : zeroe = allmin[pair[0].name].energy-point[-1]
            point.append(pair.TSstr.energy - zeroe)
            point.append(pair.TSstr.energy - zeroe)

            FS = allmin[pair[1].name].energy - zeroe
            point.append(FS)
            point.append(FS)

            tmp.append(allmin[pair[0].name])
            tmp.append(pair.TSstr)
            tmp.append(allmin[pair[1].name])

        tmp.printall('%s.arc' % (outfile))
        plotx = np.linspace(1 - 0.125,
                            len(self) + 1 + 0.125, 4 * len(self) + 2)
        return plotx, point

    def OutPair_Sorted(self, allmin, outfile="sorted"):
        for pair in self:
            pair.barrier = pair.TSstr.energy - allmin[pair[0].name].energy

        self.sort(key=lambda x: x.barrier)

        tmp = allstr()
        for pair in self:
            tmp.append(allmin[pair[0].name])
            tmp.append(pair.TSstr)
            tmp.append(allmin[pair[1].name])

        tmp.printall('%s.arc' % (outfile))
        return

    def OutPair_Sorted2(self, allmin, outfile="sorted"):
        for pair in self:
            pair.barrier = pair.TSstr.energy - allmin[pair[0].name].energy

        self.sort(key=lambda x: x.barrier)

        tmp = allstr()
        for pair in self:
            tmp.append(pair[0])
            tmp.append(pair.TSstr)
            tmp.append(pair[1])

        tmp.printall('%s.arc' % (outfile))
        return


# GibbsCorrDict = {
#     # NIST exp data
#     # https://janaf.nist.gov/
#     # - T * S + [ H(T) - H(0K) ] + kb * T * ln(P/P0)
#     "H2":
#     -500 * 145.737 * 1.036e-5 + (5.882 + 8.467) * 1.036e-2 +
#     8.61734e-5 * 500 * log(40 / 1),  # + 0.27,
#     "H2O":
#     -500 * 206.534 * 1.036e-5 + (6.925 + 9.904) * 1.036e-2 +
#     8.61734e-5 * 500 * log(1 / 1),  # + 0.52,
#     "CO2":
#     -500 * 234.901 * 1.036e-5 + (8.305 + 9.364) * 1.036e-2 +
#     8.61734e-5 * 500 * log(10 / 1),  # + 0.31,
#     "CO":
#     -500 * 212.831 * 1.036e-5 + (5.931 + 8.671) * 1.036e-2 +
#     8.61734e-5 * 500 * log(10 / 1) - 0.44,  # + 0.13,
#     "HCHO":
#     -500 * 239.102 * 1.036e-5 + (7.841 + 10.021) * 1.036e-2 +
#     8.61734e-5 * 500 * log(1 / 1),  # + 0.70,
#     # Gibbs Correction calculated by vaspkit. Remove ZPE, which is calculated by hNNCalPhono. 273.15K, 1 atm.
#     "HCOOH":
#     -0.320143 - 0.886839 + 8.61734e-5 * 500 * log(1 / 1),
#     "CH3OH":
#     0.202312 - 1.355580 + 8.61734e-5 * 500 * log(1 / 1),
# }
# GibbsCorrDict["HCOOH"] = 0
# GibbsCorrDict["HCHO"] = 0
# GibbsCorrDict["H2-HadsHads"] = GibbsCorrDict["H2"] + -0.5

# NameList = ["CO2", "CO", "H2O", "H2", "HCHO", "CH3OH", "HCOOH"]


class allstr(allstr_ecfp):
    # def AddGibbsEnergy(self, LGibbs, strType="allstr", filemode=1):
    #     """"""
    #     print("Fake Energy: Gibbs added. Temperature: %s" % LGibbs)
    #     # print("GibbsCorrection")
    #     print(GibbsCorrDict)
    #     if LGibbs <= 0:
    #         for struc in self:
    #             struc.trueEnergy = struc.energy
    #     else:
    #         """Calculate Gas Entropy by experimental data"""
    #         for struc in self:
    #             struc.trueEnergy = struc.energy

    def CountGasSpices(self, sminame):
        fragments = sminame.split('.')
        spicesDict = {}
        for frag in fragments:
            # if frag in NameList:
            if frag in GibbsCorrDict.keys():

                formula = frag
                spicesDict[
                    formula] = 1 if formula not in spicesDict else spicesDict[
                        formula] + 1
        return spicesDict

    def CountAllGasSpices(self, strType="allstr", filemode=1):
        if strType == "allstr":
            for struc in self:
                struc.gasMole = self.CountGasSpices(struc.sminame)
        if strType == "pair":
            if filemode == 1:
                for iStr in range(0, len(self), 3):
                    nIS = self.CountGasSpices(self[iStr].sminame)
                    nFS = self.CountGasSpices(self[iStr + 1].sminame)
                    nTS = {}
                    for name in nIS.keys():
                        if name in nFS.keys():
                            nTS[name] = min(nIS[name], nFS[name])

                    self[iStr].gasMole = nIS
                    self[iStr + 2].gasMole = nTS
                    self[iStr + 1].gasMole = nFS

            if filemode == 2:
                for iStr in range(0, len(self), 3):
                    nIS = self.CountGasSpices(self[iStr].sminame)
                    nFS = self.CountGasSpices(self[iStr + 2].sminame)
                    nTS = {}
                    for name in nIS.keys():
                        if name in nFS.keys():
                            nTS[name] = min(nIS[name], nFS[name])

                    self[iStr].gasMole = nIS
                    self[iStr + 1].gasMole = nTS
                    self[iStr + 2].gasMole = nFS


class AllPair(Pair):
    def __init__(self, LConly=0):
        list.__init__(self)
        self.allname = {}
        self.allnameid = {}
        self.allstrbyname = {}
        self.allstrmin = {}
        self.LGibbs = -1
        self.addmin = allstr()
        self.rpdict = {}
        self.multi = False
        self.multiSmiNameDict = {}  # {normalized sminame: normalized ecfpname}
        self.multiEcfpNameDict = {}  # {ecfpname: normalized ecfpname}
        # LGibbs == -1 mean Gibbs Correct Done
        # LGibbs == 0  mean No Gibbs Correct need be done(No calculation of zpe)
        self.LConly = LConly

    def readfile(self,
                 filename='allgoodsect.arc',
                 filemode=1,
                 Lallmininput=0,
                 LGibbs=-1,
                 allminfile='allname.arc',
                 pairsortflag=True):
        _tmp = allstr()

        _tmp.readfile(filename)
        print("read allgoodsect.arc done")

        _tmp.GetAllSminame_FromEcfp(numproc=24, colorflag=0)
        print("Gen smi  done")
        # _tmp.GetTmpFakebmx()
        # print("Gen fake BondMatrix done")
        _tmp.GetAllECFPname(numproc=24, LConly=self.LConly)
        print("Gen  ECFP done")
        _tmp.CountAllGasSpices(strType="pair", filemode=filemode)
        self.LGibbs = LGibbs
        # _tmp.AddGibbsEnergy(LGibbs, strType="pair", filemode=filemode)
        # print("Gibbs correction done")

        for i in range(0, len(_tmp), 3):
            if filemode == 1:
                # if (not _tmp[i].strflag) or (not _tmp[i + 1].strflag):
                #     continue
                # if (not _tmp[i].sminame) or (not _tmp[i + 1].sminame):
                #     continue
                # _tmp[i].GenECFPname()
                # _tmp[i+1].GenECFPname()
                _tmp[i].name = _tmp[i].ECFPname
                _tmp[i + 1].name = _tmp[i + 1].ECFPname
                TS = _tmp[i + 2].energy
                self.append(
                    Pair(_tmp[i], _tmp[i + 1], TS, _tmp[i + 2], pairsortflag))
            elif filemode == 2:
                # if (not _tmp[i].Lminstr) or (not _tmp[i+2].Lminstr): continue
                # if (not _tmp[i].sminame) or (not _tmp[i + 2].sminame):
                #     continue
                # if (not _tmp[i].strflag) or (not _tmp[i + 2].strflag):
                #     continue
                # _tmp[i].GenECFPname()
                # _tmp[i+2].GenECFPname()
                _tmp[i].name = _tmp[i].ECFPname
                _tmp[i + 2].name = _tmp[i + 2].ECFPname

                TS = _tmp[i + 1].energy
                self.append(
                    Pair(_tmp[i], _tmp[i + 2], TS, _tmp[i + 1], pairsortflag))

        self.npair = len(self)

        self.Lallmin = Lallmininput

        self.allname = {}
        self.allnameid = {}
        self.allstrbyname = {}
        _allmin = allstr()
        if Lallmininput == 1:

            _allmin.readfile(allminfile)

            _allmin.GetAllSminame_FromEcfp(numproc=24, colorflag=0)
            # _allmin.GetTmpFakebmx()
            _allmin.GetAllECFPname(LConly=self.LConly)
            _allmin.CountAllGasSpices(strType="allstr")
            # _allmin.AddGibbsEnergy(LGibbs, strType="allstr")

        if Lallmininput == 2:
            allminfile = 'addmin.arc'
            self.addallmin = allminfile

            for pair in self:
                for Str in pair:
                    if Str.name not in self.allname.keys():
                        nameid = len(self.allname)
                        self.allname[Str.name] = nameid
                        self.allnameid[nameid] = Str.name
                        self.allstrbyname[Str.name] = [Str]
                    else:
                        self.allstrbyname[Str.name].append(Str)

            _allmin.readfile(self.addallmin)
            _allmin.GetAllSminame_FromEcfp(numproc=24, colorflag=0)
            _allmin.GetAllECFPname(numproc=24, LConly=self.LConly)
            _allmin.CountAllGasSpices(strType="allstr")
            self.addmin = _allmin
            # _allmin.AddGibbsEnergy(LGibbs, strType="allstr")
            self.Lallmin = 1

        for Str in _allmin:
            Str.name = Str.ECFPname
            if Str.name not in self.allname.keys():
                nameid = len(self.allname)
                self.allname[Str.name] = nameid
                self.allnameid[nameid] = Str.name
                self.allstrbyname[Str.name] = [Str]
            else:
                self.allstrbyname[Str.name].append(Str)

        self.allstrmin = {}
        for name in self.allname.keys():
            self.allstrmin[name] = min(self.allstrbyname[name],
                                       key=lambda x: x.energy)
        self.GenAllStrMin()
        # self.allstrminTrueEnergy = {}
        # for name in self.allname.keys():
        #     self.allstrminTrueEnergy[name] = min(self.allstrbyname[name],
        #                                          key=lambda x: x.trueEnergy)

        print("all read done")

    def GibbsCorrectAll(self):
        for pair in self:
            pair[0].trueEnergy = pair[0].energy
            pair.TSstr.trueEnergy = pair.TSstr.energy
            pair[1].trueEnergy = pair[1].energy
        for struc in self.allstrmin.values():
            struc.trueEnergy = struc.energy
        for struc in self.addmin:
            struc.trueEnergy = struc.energy

        if self.LGibbs < 0:
            return

        elif self.LGibbs >= 0:
            temp = self.LGibbs
            print(
                "Gibbs Correct On LGibbs=%s: NN-ZPE + exp Entropy (only have 500 oC now)"
                % temp)
            print(GibbsCorrDict)
            print(
                "Energy of Structure without calculating ZPE will be set to be Zero"
            )
            for pair in self:
                pair[0].energy = self.GibbsCorrSingle(pair[0], temp, type=0)
                pair.TSstr.energy = self.GibbsCorrSingle(pair.TSstr,
                                                         temp,
                                                         type=1)
                pair[1].energy = self.GibbsCorrSingle(pair[1], temp, type=0)
                pair.TS = pair.TSstr.energy
            for struc in self.allstrmin.values():
                struc.energy = self.GibbsCorrSingle(struc, temp, type=0)
            for struc in self.addmin:
                struc.energy = self.GibbsCorrSingle(struc, temp, type=0)

    def GibbsCorrSingle(self, struc, temp, type=0, zpe=None):
        if zpe is not None and not hasattr(struc, "zpe"):
            # print("Struc Has no ZPE: Energy " + str(struc.energy))
            return 0
        if zpe is None:
            struc.zpe = 0

        if temp == 0:
            return struc.trueEnergy + struc.zpe

        energy = struc.trueEnergy
        for k, v in struc.gasMole.items():
            if k in GibbsCorrDict:
                energy += GibbsCorrDict[k] * v
        return energy

    # def Correct2Ref(self):
    #     for pair in self:
    #         pair[0].trueEnergy = pair[0].energy
    #         pair.TSstr.trueEnergy = pair.TSstr.energy
    #         pair[1].trueEnergy = pair[1].energy
    #     for struc in self.allstrmin.values():
    #         struc.trueEnergy = struc.energy
    #     for struc in self.addmin:
    #         struc.trueEnergy = struc.energy

    #     for pair in self:
    #         pair[0].energy = self.Correct2RefSingle(pair[0])
    #         pair.TSstr.energy = self.Correct2RefSingle(pair.TSstr)
    #         pair[1].energy = self.Correct2RefSingle(pair[1])
    #         pair.TS = pair.TSstr.energy
    #     for struc in self.allstrmin.values():
    #         struc.energy = self.Correct2RefSingle(struc)
    #     for struc in self.addmin:
    #         struc.energy = self.Correct2RefSingle(struc)

    # def Correct2RefSingle(self, struc):
    #     H2 = -6.774287 - 0.606
    #     CO2 = -22.955 - 1.034
    #     H2O = -14.219 - 0.895

    #     struc.calAtomnum()
    #     content = dict(zip(struc.elenameList, struc.natompe))
    #     # print(content)
    #     energy = struc.trueEnergy
    #     energy -= CO2 * content["C"]
    #     content["C"] -= 0
    #     content["O"] -= content["C"] * 2

    #     energy -= H2O * content["O"]
    #     content["O"] -= 0
    #     content["H"] -= content["O"] * 2

    #     energy -= H2 / 2 * content["H"]
    #     return energy

    def GenAllStrMin(self):
        self.allstrmin = {}
        for name in self.allname.keys():
            self.allstrmin[name] = min(self.allstrbyname[name],
                                       key=lambda x: x.energy)

    def JoinAllPair(self, allAnalyzer):
        if self.LGibbs < 0:
            H2 = -6.774276  # - 500 * 145.737 * 1.036e-5 + 1.3806e-23 * 500 * log( 33.273 / 1)
            H2O = -14.219145  # - 500 * 206.534 * 1.036e-5 + 1.3806e-23 * 500 * log(0.069 / 1)
        else:
            H2 = -6.774276 + GibbsCorrDict["H2"]
            H2O = -14.219145 + GibbsCorrDict["H2O"]


# {'H2': -0.4473203776244291, 'H2O': -0.8954976800000001, 'CO2': -0.9345255468736936,
#  'CO': -0.8519770668736936, 'HCHO': -1.05349804, 'HCOOH': -1.206982, 'CH3OH': -1.153268, 'H2-HadsHads': -0.34732037762442913}

        BallanceCorrection = [
            H2 * 2,
            H2O,
            H2,
        ]
        # BallanceCorrection = [
        #     0,
        #     0,
        #     0,
        # ]
        BallanceSminame = [".H2.H2", ".H2O", ".H2"]
        # BallanceSminame = ["", "", ""]
        "H2-HadsHads & H2 are thought  as same"

        # for later print: collect the connected node Structure
        jointStruc = allstr()
        jointWarnStruc = allstr()

        for index, analyzer in enumerate(allAnalyzer):
            print(len(analyzer))
            # need : allstrmin/rpdict{ecfpname:pairll}/allname{ecfpname:id}
            for struc in analyzer.allstrmin.values():
                """min struc processing"""

                fragments = [
                    re.sub(r"-.*ads", "", frag)
                    for frag in struc.sminame.split('.')
                ]

                # existFlag = False
                _normSmiName = fragments + BallanceSminame[index].split(".")
                # _normSmiName = [
                #     _n if _n != "H2" else "H.H" for _n in _normSmiName
                # ]
                normSmiName = ".".join(sorted(_normSmiName))

                if normSmiName not in self.multiSmiNameDict:
                    # print(normSmiName)
                    """IF the structure is new to the joint Graph/Net """
                    nameid = len(self.allname)
                    self.allname[struc.name] = nameid
                    self.allnameid[nameid] = struc.name

                    # self.multiSmiNameDict = {}  # {normalized sminame: normalized ecfpname}
                    # self.multiEcfpNameDict = {}  # {ecfpname: normalized ecfpname}

                    self.multiSmiNameDict[normSmiName] = struc.name
                    self.multiEcfpNameDict[struc.name] = struc.name
                    self.allstrmin[struc.name] = struc
                    if not hasattr(struc, "ballanced"):
                        struc.energy = struc.energy + BallanceCorrection[index]
                        struc.sminame = struc.sminame + BallanceSminame[index]
                        struc.ballanced = True
                else:
                    """
                    IF the structure EXISTED in the joint Graph/Net
                    join them / add it to  multiEcfpNameDict with {ecfpname: normalized ecfpname}
                    """
                    self.multiEcfpNameDict[
                        struc.name] = self.multiSmiNameDict[normSmiName]
                    self.allname[struc.name] = self.allname[
                        self.multiSmiNameDict[normSmiName]]
                    if not hasattr(struc, "ballanced"):
                        struc.energy = struc.energy + BallanceCorrection[index]
                        struc.sminame = struc.sminame + BallanceSminame[index]
                        struc.ballanced = True

                    jointStruc.append(
                        self.allstrmin[self.multiSmiNameDict[normSmiName]])
                    jointStruc.append(struc)
                    """EXAMIN if the energy difference of the joint Node is too large"""
                    if abs(struc.energy - self.allstrmin[
                            self.multiSmiNameDict[normSmiName]].energy) > 0.1:
                        print(
                            "Warning: Joint Node difference > 0.1: %.3f\t%.3f\t: %s --- %s "
                            %
                            (struc.energy, self.allstrmin[
                                self.multiSmiNameDict[normSmiName]].energy,
                             struc.sminame, self.allstrmin[
                                 self.multiSmiNameDict[normSmiName]].sminame))
                        jointWarnStruc.append(
                            self.allstrmin[self.multiSmiNameDict[normSmiName]])
                        jointWarnStruc.append(struc)
            """Ballance pair energy"""
            for key, value in analyzer.rpdict.items():
                for pair in value:
                    pair[0].name = self.multiEcfpNameDict[pair[0].name]
                    pair[1].name = self.multiEcfpNameDict[pair[1].name]
                    if not hasattr(pair[0], "ballanced"):
                        pair[0].energy += BallanceCorrection[index]
                        pair[0].sminame += BallanceSminame[index]
                        pair[0].ballanced = True
                    if not hasattr(pair[1], "ballanced"):
                        pair[1].energy += BallanceCorrection[index]
                        pair[1].sminame += BallanceSminame[index]
                        pair[1].ballanced = True
                    if not hasattr(pair.TSstr, "ballanced"):
                        pair.TSstr.energy += BallanceCorrection[index]
                        pair.TSstr.sminame += BallanceSminame[index]
                        pair.TSstr.ballanced = True
                    pair.TS = pair.TSstr.energy

            # for key, value in analyzer.rpdict.items():
            #     print("%s %s %s" % (key, self.multiEcfpNameDict[key],
            #                         analyzer.allstrmin[key].sminame))
            #     key = self.multiEcfpNameDict[key]
            #     for pair in value:
            #         tar = pair.GetFS(key).name
            #         sminame = pair.GetFS(key).sminame
            #         print("   %s %s" % (tar, sminame))

            for key, value in analyzer.rpdict.items():
                targetIS = self.multiEcfpNameDict[key]
                if targetIS not in self.rpdict:
                    """when IS is newly added"""
                    self.rpdict[targetIS] = value
                else:
                    """when IS is existed """
                    existFS = [
                        pair.GetFS(targetIS).name
                        for pair in self.rpdict[targetIS]
                    ]

                    for pair in value:
                        targetFS = self.multiEcfpNameDict[pair.GetFS(
                            targetIS).name]
                        if targetFS not in existFS:
                            self.rpdict[targetIS].append(pair)

            # print(len(analyzer.rpdict.keys()))
        # for i, struc in enumerate(jointStruc):
            # print(i + 1, struc.energy, struc.trueEnergy, struc.sminame)
        jointStruc.printall("jointNode.arc")
        jointWarnStruc.printall("jointWarn.arc")
        pass
        self.Lallmin = 1

    def ReadPath(self,
                 filename,
                 filemode=2,
                 Lallmininput=0,
                 allminfile='allname.arc'):
        self.readfile(filename,
                      filemode=filemode,
                      Lallmininput=Lallmininput,
                      allminfile=allminfile,
                      pairsortflag=False)
        outpath = Path()
        for pair in self:
            outpath.append(pair)
        return outpath

    def AllName(self):
        self.allname = {}
        self.allnameid = {}
        self.allstrbyname = {}
        for pair in self:
            for Str in pair:
                if Str.name not in self.allname.keys():
                    nameid = len(self.allname)
                    self.allname[Str.name] = nameid
                    self.allnameid[nameid] = Str.name
                    self.allstrbyname[Str.name] = [Str]
                else:
                    self.allstrbyname[Str.name].append(Str)

        if self.Lallmin == 2:
            _allmin = allstr()
            _allmin.readfile(self.addallmin)
            _allmin.GetAllSminame_FromEcfp(numproc=8, colorflag=0)
            _allmin.GetAllECFPname(LConly=self.LConly)
            _allmin.AddGibbsEnergy(self.LGibbs, strType="allstr")

            for Str in _allmin:
                Str.name = Str.ECFPname
                if Str.name not in self.allname.keys():
                    nameid = len(self.allname)
                    self.allname[Str.name] = nameid
                    self.allnameid[nameid] = Str.name
                    self.allstrbyname[Str.name] = [Str]
                else:
                    self.allstrbyname[Str.name].append(Str)

        self.allstrmin = {}
        for name in self.allname.keys():
            self.allstrmin[name] = min(self.allstrbyname[name],
                                       key=lambda x: x.energy)
        # self.allstrminTrueEnergy = {}
        # for name in self.allname.keys():
        #     self.allstrminTrueEnergy[name] = min(self.allstrbyname[name],
        #                                          key=lambda x: x.trueEnergy)

        self.Lallmin = 1

    def GenPairDict(self):
        self.allpair = {}
        for pair in self:
            if pair.name not in self.allpair.keys():
                self.allpair[pair.name] = []
                self.allpair[pair.name].append(pair)
            else:
                self.allpair[pair.name].append(pair)

        self.minpair = {}
        for name in self.allpair.keys():
            self.minpair[name] = min(self.allpair[name], key=lambda x: x.TS)

    def PrintPair(self, name1, name2, printTSstr=True):
        if self.Lallmin != 1:
            self.AllName()
        self.GenPairDict()
        namelist = [name1, name2]
        if self.multi:
            namelist = [
                self.multiEcfpNameDict[name1], self.multiEcfpNameDict[name2]
            ]
        namelist.sort()
        targetname = '%s' % (namelist[0] + namelist[1])

        print(self.allpair.keys())
        self.allpair[targetname].sort(key=lambda x: x.TS)
        _tmp = allstr()
        for pair in self.allpair[targetname]:
            # print (pair.TS)
            _tmp.append(pair[0])
            if printTSstr:
                _tmp.append(pair.TSstr)
            _tmp.append(pair[1])
        _tmp.printall('selectpair.arc')

    def PrintPair_list(self, filename, printTSstr=True):
        if self.Lallmin != 1:
            self.AllName()

        self.GenPairDict()

        _tmp = allstr()
        _tmp.readfile(filename)
        _tmp.GetAllSminame_FromEcfp(colorflag=0)
        _tmp.GetAllECFPname(LConly=self.LConly)

        out = allstr()
        for i in range(0, len(_tmp), 2):
            name1 = _tmp[i].ECFPname
            name2 = _tmp[i + 1].ECFPname

            namelist = [name1, name2]
            namelist.sort()
            targetname = '%s' % (namelist[0] + namelist[1])

            try:
                self.allpair[targetname].sort(key=lambda x: x.TS)
                for pair in self.allpair[targetname]:
                    # print (pair.TS)
                    out.append(pair[0])
                    if printTSstr:
                        out.append(pair.TSstr)
                    out.append(pair[1])
            except Exception:
                print('no pair %d' % (i / 2))

        out.printall('selectpair_list.arc')

    def GenPairdict_byname(self):
        if self.Lallmin != 1:
            self.AllName()

        self.GenPairDict()
        self.rpdict = {}
        for name in self.allname.keys():
            self.rpdict[name] = []
            for pair in self.minpair.values():
                if pair[0].name == name or pair[1].name == name:
                    self.rpdict[name].append(pair)
            self.rpdict[name].sort(key=lambda X: X.TS)
        # print self.rpdict
        pass

    def OutPairdict_byname(self, name, printTSstr=True):
        # print (name)
        _tmp = allstr()
        for pair in self.rpdict[name]:
            if pair[0].name == name:
                # print (self.allname[pair[1].name],pair[1].name,pair.TS)
                _tmp.append(pair[0])
                if printTSstr:
                    _tmp.append(pair.TSstr)
                _tmp.append(pair[1])
            else:
                # print (self.allname[pair[0].name],pair[0].name,pair.TS)
                _tmp.append(pair[1])
                if printTSstr:
                    _tmp.append(pair.TSstr)
                _tmp.append(pair[0])
        _tmp.printall('spepair.arc')

    def Outallpair(self, pts=1):
        _tmp = allstr()
        # csv = ["ISname, FSname, TSEnergy, TSZpe, TSGibbs, TSGasMole\n"]
        for pair in self.minpair.values():
            # csv.append(", ".join([
            #     pair[0].sminame, pair[1].sminame,
            #     str(pair.TSstr.trueEnergy),
            #     str(pair.TSstr.zpe),
            #     str(pair.TSstr.energy),
            #     str(pair.TSstr.gasMole)
            # ]) + "\n")
            _tmp.append(pair[0])
            if pts:
                _tmp.append(pair.TSstr)
            _tmp.append(pair[1])
        _tmp.printall('allpair.arc')
        # with open("allpair.csv", "w") as fp:
        #     fp.writelines(csv)

    def connect_path_all(self, name1, ppath, depth, maxdepth=10):
        if depth > 9 or len(name1) == 0:
            return

        pmid = []
        pmidpath = []
        for i, pair in enumerate(self.rpdict[name1]):
            # if pair[0].name == name2 or pair[1].name== name2:
            if i < 10:
                # if pair[0].name == name1 and (pair[1].name not in ppath.namechain ) and (pair.TS - pair[0].energy) < 2.5:
                if (name1 in pair[0].name) and (
                        pair[1].name not in ppath.namechain
                ) and (pair.TS - pair[0].energy) < 2.5 and (
                        pair.TS - self.zeroenergy) < 3:
                    pmid.append(pair[1].name)
                    # print (pair.TS -self.zeroenergy)
                    _ppath = Path()
                    _ppath.copy(ppath)
                    _ppath.append(pair)
                    _ppath.namechain.append(pair[1].name)
                    pmidpath.append(_ppath)
                    self.searchtrace.append(pair[1].name)
                    name2 = pair[1].name
                    if name2 not in self.AllPath_fromspeStr.keys():
                        self.AllPath_fromspeStr[name2] = [_ppath]
                    else:
                        self.AllPath_fromspeStr[name2].append(_ppath)

                # elif pair[1].name == name1 and (pair[0].name not in ppath.namechain ) and (pair.TS - pair[1].energy) < 2.5:
                elif (name1 in pair[1].name) and (
                        pair[0].name not in ppath.namechain
                ) and (pair.TS - pair[1].energy) < 2.5 and (
                        pair.TS - self.zeroenergy) < 3:
                    pmid.append(pair[0].name)
                    # print (pair.TS -self.zeroenergy)
                    _ppath = Path()
                    _ppath.copy(ppath)
                    _ppath.append(pair)
                    _ppath.namechain.append(pair[0].name)
                    pmidpath.append(_ppath)
                    self.searchtrace.append(pair[0].name)
                    name2 = pair[0].name
                    if name2 not in self.AllPath_fromspeStr.keys():
                        self.AllPath_fromspeStr[name2] = [_ppath]
                    else:
                        self.AllPath_fromspeStr[name2].append(_ppath)

        # print 'cycle',depth
        # for name in pmid:
        #    print name
        depth = depth + 1
        if len(pmid) > 0:
            pdepth = [depth] * len(pmid)
            return map(self.connect_path_all, pmid, pmidpath, pdepth)
        else:
            return

    def connect_path(self, name1, name2, ppath, depth):
        if depth > 15 or len(name1) == 0:
            return

        # print(ppath)
        #        if len(name1)!= 1:
        #            r_name1=[name1[0]]
        #            name1.pop(0)
        #            return self.connect_path(r_name1,name2,ppath,depth) ,self.connect_path(name1,name2,ppath,depth)
        #        name1 =name1[0]
        pmid = []
        pmidpath = []
        for i, pair in enumerate(self.rpdict[name1]):
            # if pair[0].name == name2 or pair[1].name== name2:
            if name2 in pair[0].name or name2 in pair[1].name:
                _ppath = Path()
                _ppath.copy(ppath)
                _ppath.append(pair)
                _ppath.namechain.append(pair[1].name)
                self.allpath.append(_ppath)
            elif i < 10:
                # if pair[0].name == name1 and (pair[1].name not in ppath.namechain ) and (pair.TS - pair[0].energy) < 2.5:
                if ((name1 in pair[0].name) and
                    (pair[1].name not in ppath.namechain)) and (
                        (pair.TS - pair[0].energy) < 2.5 and
                        (pair.TS - self.zeroenergy) < 2.5):
                    # print (pair.TS -self.zeroenergy)
                    pmid.append(pair[1].name)
                    _ppath = Path()
                    _ppath.copy(ppath)
                    _ppath.append(pair)
                    _ppath.namechain.append(pair[1].name)
                    pmidpath.append(_ppath)
                    self.searchtrace.append(pair[1].name)

                # elif pair[1].name == name1 and (pair[0].name not in ppath.namechain ) and (pair.TS - pair[1].energy) < 2.5:
                elif (name1 in pair[1].name) and (
                        pair[0].name not in ppath.namechain
                ) and (pair.TS - pair[1].energy) < 2.5 and (
                        pair.TS - self.zeroenergy) < 2.5:
                    # print (pair.TS -self.zeroenergy)
                    pmid.append(pair[0].name)
                    _ppath = Path()
                    _ppath.copy(ppath)
                    _ppath.append(pair)
                    _ppath.namechain.append(pair[0].name)
                    pmidpath.append(_ppath)
                    self.searchtrace.append(pair[0].name)
        # print 'cycle',depth
        # for name in pmid:
        #    print name
        depth = depth + 1
        if len(pmid) > 0:
            pname2 = [name2] * len(pmid)
            pdepth = [depth] * len(pmid)
            for i in range(0, len(pmid)):
                self.connect_path(pmid[i], pname2[i], pmidpath[i], pdepth[i])

            # return list(map(self.connect_path, pmid, pname2, pmidpath, pdepth))
            # return self.connect_path(pmid, pname2, pmidpath, pdepth)
        else:
            return

    def GetTarget(self, file, substrate=None):
        _tmp = allstr()
        _tmp.readfile(file)

        if substrate:
            for struc in _tmp:
                struc.atom = [at for at in struc.atom if at.ele < 20]
                struc.atom.extend(substrate[0].atom)
                struc.sortatombyele()
                struc.calAtomnum()

        _tmp.GetAllSminame_FromEcfp(colorflag=0)
        # _tmp.GenCanoicalSmiles()
        for struc in _tmp:
            struc.GenECFPname(LConly=self.LConly)
            print(struc.ECFPname)

        name1 = _tmp[0].ECFPname
        name2 = _tmp[1].ECFPname

        return name1, name2

    def GetSingleTarget(self, file):
        _tmp = allstr()
        _tmp.readfile(file)
        _tmp.GetAllSminame_FromEcfp()
        _tmp[0].GenECFPname(LConly=self.LConly)
        # _tmp.GenCanoicalSmiles()
        name1 = _tmp[0].ECFPname
        return name1

    def OutInfo(self):
        f = open('allname', 'w')
        f.write('allname appear in pair\n')
        f.write('  Id   Energy      Npair    Name \n')
        _tmp = allstr()
        for i in range(len(self.allname)):
            f.write('%4d  %12.6f  %4d %-50s\n' %
                    (i, self.allstrmin[self.allnameid[i]].energy,
                     len(self.rpdict[self.allnameid[i]]), self.allnameid[i]))
            _tmp.append(self.allstrmin[self.allnameid[i]])
        _tmp.printall('allname.arc')

        f.close()

    def OutInfo_script(self):
        self.OutInfo()
        os.system('cat allname')
        # print ('  Id  Energy   Name      Npairlink')
        # for i in range(len(self.allname)):
        #    print ('%4d    %50s     %d'%(i,self.allnameid[i], len(self.rpdict[self.allnameid[i]])))

    def GetAllPathfromSpestr(self, name1):
        barrierdict = {}
        self.lowpathall = []
        for name in self.allname.keys():
            if name != name1:
                self.findpath(name1, name, output=False)
                if len(self.allpath) == 0:
                    continue
                barrierdict[name] = self.allpath[0][-1][1]
                self.lowpathall.append(
                    [name, self.allpath[0].barrier, self.allpath[0].score])
        self.lowpathall.sort(key=lambda X: X[2])
        return self.lowpathall

    def OutAllPathfromSpestr(self, name1):
        self.GetAllPathfromSpestr(name1)
        # self.findpath_all(name1)
        _tmp = allstr()
        for item in self.lowpathall:
            name = item[0]
            print(item[0], item[1])
            _tmp.append(self.allstrmin[name])
        _tmp.printall('alllow.arc')

    def findpath_all(self, name1, output=True, printstring=False):
        self.zeroenergy = self.allstrmin[name1].energy
        self.AllPath_fromspeStr = {}
        self.searchtrace = [name1]
        ppath = Path()
        ppath.namechain = [name1]
        depth = 0

        barrierdict = {}
        self.lowpathall = []
        self.connect_path_all(name1, ppath, depth)
        for name in self.AllPath_fromspeStr.keys():
            allpath = self.AllPath_fromspeStr[name]
            for path in allpath:
                path.EndtoendSort(name1)
                path.GetBarrier_list(self.allstrmin)
            allpath.sort(key=lambda x: x.score)
            barrierdict[name] = allpath[0][-1][1]
            self.lowpathall.append(
                [name, allpath[0].barrier, allpath[0].score])

        self.lowpathall.sort(key=lambda X: X[2])
        return self.lowpathall

    def findpath(self, name1, name2, output=True, printstring=False):
        self.zeroenergy = self.allstrmin[name1].energy
        self.allpath = []
        self.searchtrace = [name1]
        ppath = Path()
        ppath.namechain = [name1]
        depth = 0

        for key, value in self.rpdict.items():
            for pair in value:
                if pair[0].name == key:
                    struc1 = pair[0]
                    struc2 = pair[1]
                elif pair[1].name == key:
                    struc1 = pair[1]
                    struc2 = pair[0]
                else:
                    print("?")

            #     rmPair = []
            #     if self.AllowedReactRule(struc1.sminame,
            #                              struc2.sminame) is False:
            #         rmPair.append(pair)

            # for pair in rmPair:
            #     self.rpdict[key].remove(pair)

        self.connect_path(name1, name2, ppath, depth)
        # print self.allpath

        # print("Multiprocessing Connect Path")
        # pool = Pool(processes=numproc)
        # counter = 0
        # self.connPathQue = [(name1, name2, ppath, depth)]
        # while len(self.connPathQue) > 0:
        #     counter += 1
        #     print(counter)
        #     _tmpList1 = self.connPathQue
        #     self.connPathQue = []
        #     _tmpList2 = pool.map(self.connect_path_mp, _tmpList1)
        #     for item in _tmpList2:
        #         if item is not None:
        #             self.connPathQue += item
        # # self.connect_path_mp((name1, name2, ppath, depth))
        # pool.close()
        # pool.join()
        print("connect path done")

        # print 'start sort path'
        for path in self.allpath:
            path.EndtoendSort(name1)
            # path.GetBarrier()
            path.GetBarrier_list(self.allstrmin)
            # path.maxTS= (max(path,key=lambda x:x.TS)).TS
            # print path.maxTS

        self.allpath.sort(key=lambda x: x.score)
        # self.RultOutPath(name1)

        if output:
            self.PrintPathLink(name1, name2, printstring)
        return self.allpath

    #   inputfile = open('allpath.pkl','wb')
    #   pickle.dump(self.allpath,inputfile)
    #   inputfile.close()

    def PrintPathLink(self, name1, name2, printstring=False, flag="longScore"):
        # print 'end sort path'
        os.system('rm -rf lowestpath-*')
        f = open('lowestpath', 'w')
        f.write('lowestpath %s  %s\n' % (name1, name2))
        number = min(20, len(self.allpath))
        for i in range(number):
            path = self.allpath[i]
            path.EndtoendSort(name1)
            # path.printpath_TS('lowestpath-%d.arc'%(i+1))

            _path = allstr()
            _path.append(self.allstrmin[name1])
            for pair in path:
                _path.append(pair[0])
                _path.append(pair.TSstr)
                _path.append(pair[1])
                _path.append(self.allstrmin[pair[1].name])
            _path.printall('lowestpath-%d.arc' % (i + 1))

            if printstring:
                path.plot('lowestpath-%d.png' % (i + 1))
            if flag == "longScore":
                f.write('------lowestpath-%d--score-%13d-----\n' %
                        (i + 1, path.score))
            else:
                f.write('------lowestpath-%d--score-%13e-----\n' %
                        (i + 1, path.score))  
            for ipair, pair in enumerate(path):
                f.write(
                    '%d  IS/TS/FS  %14.8f  %14.8f  %14.8f   %s   %s\n' %
                    ((ipair + 1), pair[0].trueEnergy, pair.TSstr.trueEnergy,
                     pair[1].trueEnergy, pair[0].sminame, pair[1].sminame))
        f.close()
        #        pair[1]

        bestPath = allstr()
        for pair in self.allpath[0]:
            bestPath.extend([pair[0], pair.TSstr, pair[1]])
        bestPath.printall("bestPath.arc")

        f = open('lowestpath.csv', 'w')
        number = min(30, len(self.allpath))
        for i in range(number):
            path = self.allpath[i]
            path.EndtoendSort(name1)

            f.write("coordinate, path-%d, gibbs, %s, status\n" %
                    (i + 1, "name"))
            f.write(
                "0.0, %14.8f, %14.8f, %s\n" %
                (self.allstrmin[name1].trueEnergy,
                 self.allstrmin[name1].energy - self.allstrmin[name1].energy,
                 self.allstrmin[name1].sminame))
            for ipair, pair in enumerate(path):
                f.write(
                    '%.2f, %14.8f, %14.8f, %s, %s\n' %
                    ((ipair + 0.3, pair[0].trueEnergy, pair[0].energy -
                      self.allstrmin[name1].energy, pair[0].sminame, "IS")))
                f.write(
                    '%.2f, %14.8f, %14.8f, %s, %14.3f\n' %
                    ((ipair + 0.5, pair.TSstr.trueEnergy, pair.TSstr.energy -
                      self.allstrmin[name1].energy, pair.TSstr.sminame,
                      (pair.TSstr.trueEnergy -
                       self.allstrmin[pair[0].name].trueEnergy))))
                f.write(
                    '%.2f, %14.8f, %14.8f, %s, %s\n' %
                    ((ipair + 0.7, pair[1].trueEnergy, pair[1].energy -
                      self.allstrmin[name1].energy, pair[1].sminame, "FS")))
                f.write('%.2f, %14.8f, %14.8f, %s\n' %
                        ((ipair + 1, self.allstrmin[pair[1].name].trueEnergy,
                          self.allstrmin[pair[1].name].energy -
                          self.allstrmin[name1].energy,
                          self.allstrmin[pair[1].name].sminame)))
        f.close()

        # if self.LGibbs > 0:
        #     _tmp = allstr()
        #     # _tmp.extend(self.allstrminTrueEnergy.values())
        #     _tmp.printall('allminTrueE.arc')

    # def findpath_Dijkstra(self, name1, name2, output=True, printstring=False):
    #     self.zeroenergy = self.allstrmin[name1].energy
    #     self.allpath = []
    #     # self.searchtrace = [name1]
    #     # ppath = []
    #     # depth = 0

    #     speLen = len(self.allname.keys())
    #     print(speLen)
    #     arrayLen = ctypes.c_int(speLen)
    #     minimum4c = pointer((ctypes.c_double * speLen)(
    #         *np.array([struc.energy for struc in self.allstrmin.values()])))
    #     pairMtx4c = (ctypes.c_double * speLen * speLen)()
    #     start = ctypes.c_int(self.allname[name1])
    #     end = ctypes.c_int(self.allname[name2])
    #     returnPath = pointer((ctypes.c_int * (speLen * 1000))())
    #     returnPathLen = pointer((ctypes.c_int * 1000)())

    #     print(
    #         "{" +
    #         ", ".join([str(struc.energy)
    #                    for struc in self.allstrmin.values()]) + "}, ")
    #     for key, value in self.rpdict.items():
    #         node1Name = key
    #         for pair in value:
    #             if pair[0].name == key:
    #                 node2Name = pair[1].name
    #             elif pair[1].name == key:
    #                 node2Name = pair[0].name
    #             else:
    #                 print("?")

    #             if node1Name != node2Name:
    #                 pairMtx4c[self.allname[node1Name]][
    #                     self.allname[node2Name]] = pair.TSstr.energy
    #                 pairMtx4c[self.allname[node2Name]][
    #                     self.allname[node1Name]] = pair.TSstr.energy
    #             # print(self.allname[node1Name], self.allname[node2Name],
    #             #       pair.TS)

    #     for i in range(speLen):
    #         print("{" + ", ".join([str(e) for e in pairMtx4c[i][:]]) + "}, ")
    #     numberOfPath = ctypes.c_int(int(800))  # max to be 500
    #     pathFinder = ctypes.cdll.LoadLibrary(
    #         "/home9/shiyf/bin/kpl-path-shiyfgit/Lib_GraphPath/shiyfPathFinder.so"
    #     )
    #     # pathFinder.PathFinderBest(minimum4c, pairMtx4c, arrayLen, start, end,
    #     #                           returnPath, returnPathLen)
    #     pathFinder.PathFinderBestKth(minimum4c, pairMtx4c, arrayLen, start,
    #                                  end, numberOfPath, returnPath,
    #                                  returnPathLen)
    #     counter = 0
    #     allPathList = list(returnPath.contents)
    #     for pathLen in returnPathLen.contents:
    #         if pathLen == 0:
    #             break
    #         bpath = allPathList[counter:(counter + pathLen)]
    #         counter += pathLen
    #         path = Path()
    #         for i in range(len(bpath) - 1):
    #             n1 = self.allnameid[bpath[i]]
    #             n2 = self.allnameid[bpath[i + 1]]
    #             for pair in self.rpdict[n1]:
    #                 if pair[0].name == n2 or pair[1].name == n2:
    #                     path.append(pair)
    #                     continue
    #         self.allpath.append(path)

    #     print("TEST connect path done")

    #     for path in self.allpath:
    #         path.EndtoendSort(name1)
    #         path.GetBarrier_list(self.allstrmin)
    #     self.allpath.sort(key=lambda x: x.score)

    #     if not output:
    #         return
    #     self.PrintPathLink(name1, name2, printstring)

    def findpath_Yen(self,
                     name1,
                     name2,
                     output=True,
                     printstring=False,
                     kpath=100,
                     printFlag=False):
        self.zeroenergy = self.allstrmin[name1].energy
        # partionSum = sum([np.exp(struc.energy) for struc in self.allstrmin.values()])

        graph = nx.DiGraph()
        graph.add_nodes_from(self.allname)

        edgeList = []
        for key, value in self.rpdict.items():
            for pair in value:
                if pair[0].name == key:
                    struc1 = pair[0]
                    struc2 = pair[1]
                elif pair[1].name == key:
                    struc1 = pair[1]
                    struc2 = pair[0]
                else:
                    print("?")

                if struc1.name != struc2.name and self.AllowedReactRule(
                        struc1.sminame, struc2.sminame):
                    if pair.TSstr.energy - self.zeroenergy < 3:
                        wei = np.exp(
                            (pair.TSstr.energy - self.zeroenergy) / 0.043)
                    else:
                        continue
                    # else:
                    #     wei = np.exp(
                    #         (pair.TSstr.energy - struc1.energy) / 0.043)
                    edgeList.extend([
                        (self.allname[struc1.name], self.allname[struc2.name],
                         {
                             'weight': wei
                         }),
                    ])
        if printFlag:
            print("Yen Init Done.")
        graph.add_edges_from(edgeList)

        # pos = nx.spring_layout(graph)
        # nx.draw_networkx(graph, pos)
        # labels = nx.get_edge_attributes(graph, 'weight')
        # # nx.draw_networkx_edge_labels(graph, edge_labels=labels)
        # nx.draw_networkx_edge_labels(graph, pos, edge_labels=labels)
        # plt.savefig("graph.png")
        try:
            pathSet, costList = self.PathFinderYenAlgor(
                graph, self.allname[name1], self.allname[name2], kpath)
        except nx.NetworkXNoPath:
            return []
        self.allpath = []
        for p, c in zip(pathSet, costList):
            path = Path()
            for i in range(len(p) - 1):
                n1 = self.allnameid[p[i]]
                n2 = self.allnameid[p[i + 1]]
                for pair in self.rpdict[n1]:
                    if pair[0].name == n2 or pair[1].name == n2:
                        path.append(pair)
                        continue
            path.score = c
            self.allpath.append(path)
        if printFlag:
            print("TEST connect path done")

        for path in self.allpath:
            path.EndtoendSort(name1)

        # for path in self.allpath:
        #     path.EndtoendSort(name1)
        # path.GetBarrier_list(self.allstrmin)
        self.allpath.sort(key=lambda x: x.score)

        if output:
            self.PrintPathLink(name1, name2, printstring)
        return self.allpath

    def findpath_MK(self, name1, name2, printstring=False, output=True):
        import pathfindingMK

        print(name2, self.allstrmin[name2].sminame)
        allPath = self.findpath_Yen(name1, name2, False, False, 10)

        mkPair = list()
        for path in allPath:
            for p in path:
                if p not in mkPair:
                    mkPair.append(p)
        species, reaction = pathfindingMK.GenMKModel(self, mkPair)
                
        task = pathfindingMK.RunMkModel(species, reaction)

        if task.success:
            print(task.species)
            pathfindingMK.ScoreMKPair(name1, self, allPath,
                        task.species)
        else:
            print(
                "MKWarn: Err > 1e-9 or negatice coverage in Microkinetics")

        allPath.sort(key=lambda x: x.score)
        self.allpath = allPath
        if output:
            self.PrintPathLink(name1, name2, printstring, flag="floatScore")
        return self.allpath

    def PrintKeyPair(self, pairFile="pairSet.arc"):
        """print select pair(from lowest path of other surfaces and reaction with particular concern)"""
        self.GenMinPair(type="rpdict")
        """keyPair from pairFile"""
        keyPair = allstr()
        keyPair.readfile(pairFile)
        substrate = allstr()
        substrate.readfile("sub.arc")
        """substitute substrate of keyPair"""
        _allstr = allstr()
        _allstr.extend(keyPair)
        for struc in _allstr:
            struc.atom = [at for at in struc.atom if at.ele < 20]
            struc.atom.extend(substrate[0].atom)
            struc.sortatombyele()
            struc.calAtomnum()
        # _allstr.printall("test.arc")

        _allstr.GetAllSminame_FromEcfp(numproc=16, colorflag=0)
        print("Gen Key pair smi done")
        _allstr.GetAllECFPname(numproc=16, LConly=self.LConly)
        print("Gen key pair ECFP done")

        if self.multi is True:
            for struc in _allstr:
                struc.ECFPname = self.multiEcfpNameDict.get(
                    struc.ECFPname, struc.ECFPname)

        keyPair = AllPair()
        for i in range(0, len(_allstr), 3):
            _allstr[i].name = _allstr[i].ECFPname
            _allstr[i + 2].name = _allstr[i + 2].ECFPname
            TS = _allstr[i + 1].energy
            keyPair.append(
                Pair(_allstr[i], _allstr[i + 2], TS, _allstr[i + 1], True))
        """unique keyPair"""
        uniquePair = AllPair()
        for name in np.unique([pair.name for pair in keyPair]):
            uniquePair.append([pair for pair in keyPair
                               if pair.name == name][0])
        keyPair = uniquePair

        # _tmp = allstr()
        # for pair in keyPair:
        #     _tmp.append(pair[0])
        #     _tmp.append(pair.TSstr)
        #     _tmp.append(pair[1])
        # _tmp.printall("UniquePairSet.arc")
        """output minimum Structrue by keyPair"""
        keyName = np.unique([pair[0].name for pair in keyPair] +
                            [pair[1].name for pair in keyPair])
        keyStruc = {
            key: self.allstrmin[key]
            for key in keyName if key in self.multiEcfpNameDict
        }
        outStruc = allstr()
        outStruc.extend(list(keyStruc.values()))
        outStruc.printall("keyStruc.arc")
        limit2ads = {
            "HC(OH)2": "HC(OH)2*",
            "H3": "H3*",
        }
        for pair in keyPair:
            pair.sort(key=lambda x: x.name)
            speC1 = self.GetCompose(pair[0], limit2ads)
            speC2 = self.GetCompose(pair[1], limit2ads)
            pair.freeSite = min([speC1["*"], speC2["*"]])
            pair.reacNProd = {
                key: speC2.get(key, 0) - speC1.get(key, 0)
                for key in np.unique(list(speC1.keys()) + list(speC2.keys()))
                if speC2.get(key, 0) - speC1.get(key, 0) != 0
            }
            pair.RvReacNProd = {k: -v for k, v in pair.reacNProd.items()}

        outPair = AllPair()
        for pair in keyPair:
            if str(pair.reacNProd) in self.minpair:
                outPair.append(self.minpair[str(pair.reacNProd)])
        print(len(outPair))
        outPair.GenMinPair(type="allpair")
        outPair = list(outPair.minpair.values())
        outPair.sort(key=lambda x: str(x.reacNProd))

        _tmp = allstr()
        for pair in outPair:
            _tmp.extend([pair[0], pair.TSstr, pair[1]])
        _tmp.printall("keyPair.arc")
        for pair in outPair:
            _tmp = allstr()
            _tmp.extend([pair[0], pair.TSstr, pair[1]])
            _tmp.printall("uncm_%s_%s.arc" %
                          (pair[0].sminame, pair[1].sminame))

        outLines = []
        for i, pair in enumerate(outPair):
            # li = ["%35s, %35s, " % (pair[0].sminame, pair[1].sminame)]
            li = [str(i)]
            li.append("".join([
                ".".join([
                    "%d%s" % (-value, react)
                    for react, value in pair.reacNProd.items() if value < 0
                ]), ">>", ".".join([
                    "%d%s" % (value, react)
                    for react, value in pair.reacNProd.items() if value > 0
                ])
            ]).center(40))
            print("".join(li))
            print(pair[0].sminame, pair[1].sminame)
            # if pair.name in self.minpair:
            targetPair = pair
            li.append(
                "%.3f, %.3f, %.3f, " %
                (self.allstrmin[pair[0].name].energy, targetPair.TSstr.energy,
                 self.allstrmin[pair[1].name].energy))
            li.append("%.3f, %.3f, %.3f," %
                      (self.allstrmin[pair[0].name].trueEnergy,
                       targetPair.TSstr.trueEnergy,
                       self.allstrmin[pair[1].name].trueEnergy))
            li.append("\n")
            # else:
            #     li.append("0,         " * 6 + "\n")
            outLines.append(li)
        outLines = ["".join(li) for li in outLines]
        with open("keyInfo.csv", "w") as fp:
            fp.writelines(outLines)

    def PrintKeyPairTarget(self, filename):

        pair = allstr()
        pair.readfile(filename)
        pair.GetAllSminame_FromEcfp(colorflag=0)
        pair[0].GenECFPname(LConly=self.LConly)
        pair[1].GenECFPname(LConly=self.LConly)
        pair[0].name = pair[0].ECFPname
        pair[1].name = pair[1].ECFPname

        pair.sort(key=lambda x: x.name)
        speC1 = self.GetCompose(pair[0])
        speC2 = self.GetCompose(pair[1])
        pair.freeSite = min([speC1["*"], speC2["*"]])
        pair.reacNProd = {
            key: speC2.get(key, 0) - speC1.get(key, 0)
            for key in np.unique(list(speC1.keys()) + list(speC2.keys()))
            if speC2.get(key, 0) - speC1.get(key, 0) != 0
        }

        self.GenMinPair(type="rpdict")
        print(self.allminpair.keys())
        _tmp = allstr()
        for p in self.allminpair[str(pair.reacNProd)]:
            _tmp.extend([p[0], p.TSstr, p[1]])
        _tmp.printall("kptarget.arc")

    def GetCompose(self, struc, limit2ads={}):
        compos = list(np.unique(struc.sminame.split("."), return_counts=True))

        compos[0] = [limit2ads.get(key, key) for key in compos[0]]
        compos[0] = [re.sub(r"-.*ads", "", key) for key in compos[0]]
        compos[0] = [re.sub(r"\(", "<", key) for key in compos[0]]
        compos[0] = [re.sub(r"\)", ">", key) for key in compos[0]]
        compos[0] = [key + "*" for key in compos[0]]
        # compos = dict(zip(list(compos[0]), list(compos[1])))
        _tmp = dict()
        for index, key in enumerate(compos[0]):
            if key not in _tmp:
                _tmp[key] = compos[1][index]
            else:
                _tmp[key] += compos[1][index]
        compos = _tmp

        compos["*"] = 9 - sum(
            [key.count("*") * value for key, value in compos.items()])
        # print(compos)
        return compos

    def GetReactName(self, reacNProd):
        reac = {k: v for k, v in reacNProd.items() if v < 0}
        prod = {k: v for k, v in reacNProd.items() if v > 0}
        reac = ".".join(
            ["%d%s" % (v, k) if v != 1 else "%s" % k for k, v in reac.items()])
        prod = ".".join(
            ["%d%s" % (v, k) if v != 1 else "%s" % k for k, v in prod.items()])
        return (" >> ".join([reac, prod]), " >> ".join([reac, prod]))

    def GenMinPair(self, name1="", name2="", type="yen"):
        limit2ads = {"HC(OH)2": "HC(OH)2*", "H3": "H3*"}
        _allpair = []

        # if KSP Yen first
        if type == "yen":
            self.findpath_Yen(name1, name2, output=False, kpath=500)
            _tmp = allstr()
            for path in self.allpath[0:3]:
                for pair in path:
                    if pair not in _allpair:
                        _allpair.append(pair)
                        _tmp.append(pair[0])
                        _tmp.append(pair.TSstr)
                        _tmp.append(pair[1])
            _tmp.printall("mkm.arc")

        # for allpair
        elif type == "rpdict":
            _tmp = allstr()
            for lst in self.rpdict.values():
                for pair in lst:
                    # if pair not in _allpair and pair.TS - self.allstrmin[
                    #         name1].energy < 3:
                    if pair not in _allpair:
                        _allpair.append(pair)
                        _tmp.append(pair[0])
                        _tmp.append(pair.TSstr)
                        _tmp.append(pair[1])
            _tmp.printall("mkm.arc")

        elif type == "allpair":
            _allpair = self

        for pair in _allpair:
            pair.sort(key=lambda x: x.name)

        for struc in self.allstrmin.values():
            struc.speCount = self.GetCompose(struc, limit2ads)

        for pair in _allpair:
            if self.multi:
                speC1 = self.allstrmin[pair[0].name].speCount
                speC2 = self.allstrmin[pair[1].name].speCount
            else:
                speC1 = self.GetCompose(pair[0], limit2ads)
                speC2 = self.GetCompose(pair[1], limit2ads)
            pair.freeSite = min([speC1["*"], speC2["*"]])
            # print(speC1)
            # print(speC2)
            pair.reacNProd = {
                key: speC2.get(key, 0) - speC1.get(key, 0)
                for key in np.unique(list(speC1.keys()) + list(speC2.keys()))
                if speC2.get(key, 0) - speC1.get(key, 0) != 0
            }
            # print(pair.reacNProd)
            # print(pair[0].sminame, pair[1].sminame)

        self.minpair = {}
        for pair in _allpair:
            key = str(pair.reacNProd)
            if key not in self.minpair:
                self.minpair[key] = pair
            else:
                if pair.freeSite > self.minpair[key].freeSite:
                    self.minpair[key] = pair

        self.allminpair = {}
        for pair in _allpair:
            key = str(pair.reacNProd)
            # print(key)
            if key not in self.allminpair:
                self.allminpair[key] = [pair]
            else:
                self.allminpair[key].append(pair)

    def findpath_MKMCXX(self, name1, name2, output=True, printstring=False):
        gasPhsPressure = {
            "CO2": 0.2,
            "CO": 1.8,
            "H2": 20,
            "HCOOH": 1e-4,
            "HCHO": 1e-4,
            "CH3OH": 0.01,
            "H2O": 0.01,
            "CH4": 1e-5,
        }

        self.GenMinPair(name1, name2, type=0)

        speSet = set()
        allReacNProd = []
        for pair in self.minpair.values():
            print(pair.reacNProd)
            allReacNProd.append(pair.reacNProd)
            speSet.update(pair.reacNProd.keys())

        mkm = [
            """&settings
TYPE = SEQUENCERUN
PRESSURE = -1
DEBUG = 4
#SOLMAXSTEP = 10000
#SOLTESTFAIL = 200
#SOLCONVFAIL = 200
#SOLVERTYPE = 1
#BOOSTER = 1
MAKEPLOTS = 0
&runs
#
# Temp;	Time;	AbsTol;	RelTol
500;  2e10;	 1e-10;   1e-10
"""
        ]
        adsorbtion = """
AR; {H2*} => {*} + {H2}; 1e13; 1e13; 0;  57505
AR; {CO*} => {*} + {CO};    1e13; 1e13; 0; 9381
AR; {CO2*} => {*} + {CO2};    1e13; 1e13; 0; 98087
AR; {H2O*} => {*} + {H2O};    1e13; 1e13; 0; 73743
AR; {HCOOH*} => {*} + {HCOOH};    1e13; 1e13; 0; 101284
AR; {HCHO*} => {*} + {HCHO};    1e13; 1e13; 0; 102947
AR; {CH3OH*} => {*} + {CH3OH};    1e13; 1e13; 0; 96606
"""

        mkm.append("&compounds\n")
        for spe, pressure in gasPhsPressure.items():
            mkm.append("%s;    0;   %f\n" % (spe, pressure))
        mkm.append("*;   1;    1\n")
        for spe in speSet:
            if "*" in spe and spe != "*":
                mkm.append("%s;    1;   0\n" % (spe))

        ev2j = 96485.37
        mkm.append("&reactions\n")
        mkm.append(adsorbtion)
        for pair in self.minpair.values():
            reac = []
            prod = []
            for spe, count in pair.reacNProd.items():
                if count < 0:
                    reac.append("%d{%s}" % (-count, spe))
                if count > 0:
                    prod.append("%d{%s}" % (count, spe))
            if len(reac) == 0 or len(prod) == 0:
                continue
            mkm.append("".join([
                "AR; ", " + ".join(reac), " => ", " + ".join(prod), ";    ",
                "1e13; 1e13; %.5e; %.5e" %
                ((pair.TSstr.energy - self.allstrmin[pair[0].name].energy) *
                 ev2j,
                 (pair.TSstr.energy - self.allstrmin[pair[1].name].energy) *
                 ev2j), "\n"
            ]))

            print(
                "%.3f" %
                (pair.TSstr.energy - self.allstrmin[pair[0].name].energy),
                "%.3f" %
                (pair.TSstr.energy - self.allstrmin[pair[1].name].energy),
                "%.3f" % (pair.TSstr.trueEnergy -
                          self.allstrmin[pair[0].name].trueEnergy),
                "%.3f" % (pair.TSstr.trueEnergy -
                          self.allstrmin[pair[1].name].trueEnergy),
                pair[0].sminame, pair[1].sminame, pair.reacNProd)

        with open("CO2Hydro.mkm", "w") as fp:
            fp.writelines(mkm)
        pass

    def findpath_MKMCXX_Post(self,
                             name1,
                             name2,
                             output=True,
                             printstring=False):
        # latestRun = new_report(os.getcwd())
        # if latestRun is None:
        #     return {}
        # with open(latestRun + "/range/rates.dat") as fp:
        #     lines = fp.readlines()

        self.GenMinPair(name1, name2)

        with open("rates.dat") as fp:
            lines = fp.readlines()

        key = lines.pop(0).split("\t")
        value = lines[0].split("\t")
        absoluteRates = {}
        lines = []
        for i in range(1, len(key), 2):
            print(key[i].replace(" ", ""), value[i], value[i + 1])
            absoluteRates[key[i].replace(" ", "").replace(
                "->", "=>")] = float(value[i]) - float(value[i + 1])
            lines.append("%s    %s     %s\n" %
                         (key[i], value[i], value[i + 1]))

        print(absoluteRates)
        for pair in self.minpair.values():
            reac = []
            prod = []
            for spe, count in pair.reacNProd.items():
                if count < 0:
                    if count == -1:
                        reac.append("%s" % (spe))
                    else:
                        reac.append("%d%s" % (-count, spe))
                if count > 0:
                    if count == 1:
                        prod.append("%s" % (spe))
                    else:
                        prod.append("%d%s" % (count, spe))
            if len(reac) == 0 or len(prod) == 0:
                continue
            pair.reactString = "".join(["+".join(reac), "=>", "+".join(prod)])

            pair.weight = absoluteRates.get(pair.reactString, 0)

        weightList = sorted(self.minpair.values(),
                            key=lambda p: abs(p.weight)
                            if hasattr(p, "weight") else 0,
                            reverse=True)
        outputPair = allstr()
        outputStruc = allstr()
        maxPair = 50
        print("MaxPair: %s" % maxPair)
        for pair in weightList[0:maxPair]:
            outputPair.append(pair[0])
            outputPair.append(pair.TSstr)
            outputPair.append(pair[1])
            if self.allstrmin[pair[0].name] not in outputStruc:
                outputStruc.append(self.allstrmin[pair[0].name])
            if self.allstrmin[pair[1].name] not in outputStruc:
                outputStruc.append(self.allstrmin[pair[1].name])
            print("%.4e" % pair.weight, end="\t")
            print("  ".join([
                "%.4f" % x for x in [
                    self.allstrmin[pair[0].name].energy,
                    pair[0].energy,
                    pair.TSstr.energy,
                    pair[1].energy,
                    self.allstrmin[pair[1].name].energy,
                    self.allstrmin[pair[0].name].trueEnergy,
                    pair[0].trueEnergy,
                    pair.TSstr.trueEnergy,
                    pair[1].trueEnergy,
                    self.allstrmin[pair[1].name].trueEnergy,
                ]
            ]),
                  end="  ")
            print(pair.reactString + "\t\t" + pair[0].sminame + " => " +
                  pair[1].sminame)
        outputPair.printall("selectPair.arc")
        outputStruc.printall("selectStruc.arc")

    def findpath_MKMCXX_FromVasp(self):
        pairDir = "VASP"
        minstrcDir = "VASP_Struc"
        outputName = "CO2Hydro_VASP.mkm"
        rootdir = os.getcwd()

        os.chdir(minstrcDir)
        minstr = allstr()
        for f in os.listdir():
            minstr.readfile(f + "/best.arc")
        os.chdir(rootdir)
        minstr.GetAllSminame_FromEcfp(numproc=16, colorflag=0)
        minstr.GetAllECFPname(numproc=16, LConly=self.LConly)
        minstr = {pair.sminame: pair for pair in minstr}

        os.chdir(pairDir)
        for f in os.listdir():
            _tmpTS = allstr()
            _tmpTS.readfile(f + "/TSstr.arc")
            _tmpIF = allstr()
            _tmpIF.readfile(f + "/uncm.arc")
            _tmpIF.GetAllSminame_FromEcfp(numproc=16, colorflag=0)
            _tmpIF.GetAllECFPname(numproc=16, LConly=self.LConly)
        os.chdir(rootdir)

    def OutDotplotfile(self, name1, name2):
        f = open('circuit.gv', 'w')

        f.write('digraph G {\n')

        f.write('ranksep = 1;\n')
        f.write('ratio = auto;\n')
        # # f.write('{ 1 -> 2 -> 3 -> 26 -> 38 -> 28 -> 34 ;}\n')
        f.write('{rank = min; "%s"};\n' % self.allstrmin[name1].sminame)
        f.write('{rank = max; "%s"};\n' % self.allstrmin[name2].sminame)
        # f.write('size ="10,10";\n')
        _outPair = allstr()

        for pair in self.circuitPair:
            if pair.current > 0:
                strucIS = pair[0]
                strucFS = pair[1]
            elif pair.current < 0:
                strucIS = pair[1]
                strucFS = pair[0]
            else:
                continue
            if abs(pair.current) / self.maxCur < self.currentLimit:
                continue
            f.write('"%s" -> "%s" [penwidth = %f, weight = %f];\n' %
                    (self.allstrmin[strucIS.name].sminame,
                     self.allstrmin[strucFS.name].sminame,
                     2.0 * abs(pair.current) / self.maxCur + 0.1,
                     log(abs(pair.current) / self.maxCur)))

            _outPair.append(pair[0])
            _outPair.append(pair.TSstr)
            _outPair.append(pair[1])

        # for name in self.allname.keys():

        #     # color = int((self.allstrmin[name].energy-self.minE)/0.23)+1
        #     color = 10 - int((self.allstrmin[name].energy - self.minE) / 0.2)
        #     f.write(
        #         '"%s"  [label = " ",height= %f, width = %f,fontsize = %f, style= filled ,colorscheme = "rdylgn10", color = %d,nodesep =0.05]\n'
        #         %
        #         (self.Strname2id[name],
        #          (1 * np.sqrt(float(self.allstrlen[name]) / self.maxlen) * 3 +
        #           1),
        #          (1 * np.sqrt(float(self.allstrlen[name]) / self.maxlen) * 4.5
        #           + 1.5),
        #          (10 * np.sqrt(float(self.allstrlen[name]) / self.maxlen) +
        #           10), color))
        f.write('}\n')
        f.close()
        _outPair.printall("mainI.arc")

        graph = nx.Graph()
        graph.add_nodes_from(self.allname)

        edgeList = []

        for pair in self.circuitPair:
            wei = 1 / abs(pair.current)
            edgeList.extend([
                (pair[0].name, pair[1].name, {
                    'weight': wei
                }),
            ])
        graph.add_edges_from(edgeList)
        p = nx.dijkstra_path(graph, name1, name2)

        self.allpath = []

        path = Path()
        for i in range(len(p) - 1):
            n1 = p[i]
            n2 = p[i + 1]
            for pair in self.rpdict[n1]:
                if pair[0].name == n2 or pair[1].name == n2:
                    path.append(pair)
                    continue
        self.allpath.append(path)

        print("TEST connect path done")

        for path in self.allpath:
            path.EndtoendSort(name1)
            path.GetBarrier_list(self.allstrmin)
        self.allpath.sort(key=lambda x: x.score)

        self.PrintPathLink(name1, name2)
        # f = open('Strinfo-ECFP', 'w')
        # for id in range(len(self.Strname2id)):
        #     name = self.Strid2name[id]
        #     f.write('%10s  %10s  minE %f\n' %
        #             (self.Strname2id[name], self.allstrmin[name].sminame,
        #              self.allstrmin[name].energy))
        # f.close()

    def OutDotRsGraph(self):

        self.Strname2id = {
            name: index
            for index, name in enumerate(self.rpdict.keys())
        }
        self.Strid2name = {v: k for k, v in self.Strname2id.items()}
        # self.Strname2id = {
        #     name: index
        #     for index, name in enumerate(self.rpdict.keys())
        # }
        f = open('Strinfo-ECFP', 'w')
        for id in range(len(self.Strname2id)):
            name = self.Strid2name[id]
            f.write('%10s  %10s  minE %f\n' %
                    (self.Strname2id[name], self.allstrmin[name].sminame,
                     self.allstrmin[name].energy))
        f.close()

        self.minE = self.allstrmin[list(self.rpdict.keys())[0]].energy
        self.maxLen = max([len(self.rpdict[k]) for k in self.rpdict])

        for k in self.rpdict:
            for pair in self.rpdict[k]:
                pair.weight = 1

        # self.findpath_Yen(self.Strid2name[6],
        #                   self.Strid2name[33],
        #                   output=False)
        # for index, path in enumerate(self.allpath):
        #     for pair in path:
        #         pair.weight += 1

        self.allpair = []
        for k in self.rpdict:
            for pair in self.rpdict[k]:
                if pair not in self.allpair:
                    self.allpair.append(pair)

        f = open('DotRsGraph.gv', 'w')
        f.write('digraph G {\n')
        f.write('nodesep = 0.1;\n')
        # f.write('ratio = 0.6;\n')
        f.write("ratio = compress\nsize = \"3600,1000\"\n")
        # f.write('{ 0 -> 1 -> 5 -> 59 -> 58 -> 60 -> 23 -> 24 -> 33;}\n')
        # f.write('{ 6 -> 28 -> 27 -> 23}\n')
        f.write('{rank = min; "0" "8"};\n')
        f.write('{rank = max; "31"};\n')

        redPath = [
            "H2.CO2.H2.H2",
            "H-Hads.H-Hads.CO2.H2.H2",
            "HCOO-OadsOads.H-Hads.H2.H2",
            "HCOO-OadsOads.H-Hads.H-Hads.H-Hads.H2",
            "HCOOH.H-Hads.H-Hads.H2",
            "H2COOH-OadsOads.H-Hads.H2",
            "HCHO.OH-Oads.H-Hads.H2",
            "HCHO.H2.H2O",
            "HCHO.H-Hads.H-Hads.H2O",
            "CH3O-Oads.H-Hads.H2O",
            "CH3OH.H2O",
            "H2O.CO-Cads.H2.H2",
            "H-Hads.H-Hads.H2.CO-Cads.H2O",
            "H-Hads.H-Hads.H-Hads.H-Hads.CO-Cads.H2O",
            "H-Hads.CH=O-CadsOads.H-Hads.H-Hads.H2O",
            "HCHO.H-Hads.H-Hads.H2O",
            "CH3O-Oads.H-Hads.H2O",
            "CH3OH.H2O",
        ]
        redLine = []
        smi2Ecfp = {
            struc.sminame: struc.name
            for struc in self.allstrmin.values()
        }
        print(smi2Ecfp)
        for i in range(len(redPath) - 1):
            edge = (self.Strname2id[smi2Ecfp[redPath[i]]],
                    self.Strname2id[smi2Ecfp[redPath[i + 1]]])
            redLine.append((min(edge), max(edge)))
        # redLine = {(0, 1), (1, 6), (6, 58), (56, 58), (56, 60), (60, 62),
        #            (20, 62), (20, 21), (21, 28), (28, 32)}
        adsEnergyCut = 2.5
        tsEnergyCut = 3.0

        nodeCollect = set()
        for pair in self.allpair:
            # pairname = pair.name
            if (self.allstrmin[pair[0].name].energy > self.minE + adsEnergyCut) \
                    or (self.allstrmin[pair[1].name].energy > self.minE + adsEnergyCut):
                continue
            if (pair.TS - self.allstrmin[pair[0].name].energy > adsEnergyCut
                    and pair.TS - self.allstrmin[pair[1].name].energy >
                    adsEnergyCut) or (pair.TS - self.minE > tsEnergyCut):
                continue
            # if (self.allstrmin[pair[0].name].Del
            #         == 1) or (self.allstrmin[pair[1].name].Del == 1):
            #     continue

            if self.Strname2id[pair[0].name] == self.Strname2id[pair[1].name]:
                print('wrong pair %s %s' % (pair[0].sminame, pair[1].sminame))
                continue
            if len(self.allstrmin[pair[0].name].sminame) > 50 or len(
                    self.allstrmin[pair[1].name].sminame) > 50:
                continue
                _tmp = allstr()
                _tmp.extend(pair)
                _tmp.printall(pair[0].name + "_" + pair[0].name + ".arc")
            print(self.Strname2id[pair[0].name], self.Strname2id[pair[1].name])

            edge = (self.Strname2id[pair[0].name],
                    self.Strname2id[pair[1].name])
            print(edge)
            if (min(edge), max(edge)) in redLine:
                f.write(
                    '"%s" -> "%s" [penwidth = %f, weight = %d, arrowhead = none, arrowtail=none, color=red];\n'
                    % (
                        self.Strname2id[pair[0].name],
                        self.Strname2id[pair[1].name],
                        max(2.0, 6 * (2.5 - (pair.TSstr.energy - self.minE))),
                        # 3,
                        5))
            else:
                f.write(
                    '"%s" -> "%s" [penwidth = %f, weight = %d, arrowhead = none, arrowtail=none];\n'
                    % (
                        self.Strname2id[pair[0].name],
                        self.Strname2id[pair[1].name],
                        max(2.5, 6 * (2.5 - (pair.TSstr.energy - self.minE))),
                        # 3,
                        1))

            nodeCollect.add(self.Strname2id[pair[0].name])
            nodeCollect.add(self.Strname2id[pair[1].name])
        for name in self.allstrmin.keys():
            # if (self.allstrmin[name].Del == 1):
            #     continue

            if self.Strname2id[name] not in nodeCollect:
                continue
            if (self.allstrmin[name].energy > (self.minE + adsEnergyCut)):
                continue
            else:
                # color = int((self.allstrmin[name].energy-self.minE)/0.23)+1
                color = max(
                    10 -
                    int(max(0,
                            (self.allstrmin[name].energy - self.minE) / 0.2)),
                    1)
                fragments = [
                    re.sub(r"-.*ads", "*", frag)
                    for frag in self.allstrmin[name].sminame.split('.')
                ]
                fragments.sort()
                label = ""
                if True:
                    rowLen = 0
                    for n in fragments:
                        if rowLen < 7:
                            label += "." + n
                            rowLen += len(n) + 1
                        else:
                            label += "<br></br><br></br>" + n
                            rowLen = 0
                    label = label[1:]
                    label = label.replace(
                        "2", "<FONT POINT-SIZE='40'><SUB>2</SUB></FONT>")
                    label = label.replace(
                        "3", "<FONT POINT-SIZE='40'><SUB>3</SUB></FONT>")
                    label = "<" + label + ">"
            if self.allstrmin[name].sminame in redPath or len(
                    self.rpdict[name]) > 8 or color > 8:
                f.write(
                    '"%s"  [label = %s, height= %f, width = %f,fontsize = %f, style= filled ,colorscheme = "rdylgn10", color = %d,nodesep =0.05]\n'
                    % (self.Strname2id[name], label,
                       1 * np.sqrt(len(self.rpdict[name]) / self.maxLen) * 3,
                       1 * np.sqrt(len(self.rpdict[name]) / self.maxLen) * 4.5,
                       50, color))
            else:
                f.write(
                    '"%s"  [label = "", height= %f, width = %f,fontsize = %f, style= filled ,colorscheme = "rdylgn10", color = %d,nodesep =0.05]\n'
                    % (
                        self.Strname2id[name],
                        # " ",
                        # label,
                        1 * np.sqrt(len(self.rpdict[name]) / self.maxLen) * 3,
                        1 * np.sqrt(len(self.rpdict[name]) / self.maxLen) *
                        4.5,
                        50,
                        color))
            # f.write(
            #     '"%s"  [style= filled ,colorscheme = "rdylgn10", color = %d]\n'
            #     % (self.Strname2id[name], color))
        f.write('}\n')
        f.close()
        os.system("cat DotRsGraph.gv")

    def AllowedReactRule(self, name1, name2):
        compos1 = np.unique(name1.split("."), return_counts=True)
        compos1 = dict(zip(list(compos1[0]), list(compos1[1])))
        compos2 = np.unique(name2.split("."), return_counts=True)
        compos2 = dict(zip(list(compos2[0]), list(compos2[1])))

        # rule 1: remove  Hydrogen Combination
        if (compos2.get("H-Hads", 0) - compos1.get("H-Hads", 0) == -2 and
            (compos2.get("H2", 0) + compos2.get("H2-HadsHads", 0) -
             compos1.get("H2", 0) - compos1.get("H2-HadsHads", 0) == 1)):
            return False

        # # rule 2: remove Hydrogen Eley-Ridel Addition Reaction
        # if (compos2.get("H-Hads", 0) - compos1.get("H-Hads", 0) == 1 and
        #     (compos2.get("H2", 0) + compos2.get("H2-HadsHads", 0) -
        #      compos1.get("H2", 0) - compos1.get("H2-HadsHads", 0) == -1)):
        #     return False

        # # # rule 3: remove CO contain route:
        # # if (compos2.get("CO-Cads", 0) == 1 or compos2.get("CO-Oads", 0) == 1
        # #         or compos1.get("CO-Cads", 0) == 1
        # #         or compos1.get("CO-Oads", 0) == 1):
        # #     return False

        # rule 4: Hydrogen or H-Hads are not allowed to be Gernerated alone:
        # # inlude rule 1:
        # if (compos2.get("H-Hads", 0) - compos1.get("H-Hads", 0) == 1 or
        #     (compos2.get("H2", 0) + compos2.get("H2-HadsHads", 0) -
        #      compos1.get("H2", 0) - compos1.get("H2-HadsHads", 0) == 1)):
        #     return False

        # # remove too high H coverage
        # if compos2.get("H-Hads", 0) > 2 or compos1.get("H-Hads", 0) > 2:
        #     return False

        # rule for formate: must contain formate species
        # if not ((compos1.get("HCOO", 0) or compos1.get("HCOO-Oads", 0)
        #          or compos1.get("HCOO-OadsOads", 0)) and
        #         (compos2.get("HCOO", 0) or compos2.get("HCOO-Oads", 0)
        #          or compos2.get("HCOO-OadsOads", 0))):
        #     return False

        # rule for path Figure:
        # if compos1.get("CH=O-Cads", 0) or compos2.get("CH=O-Cads", 0):
        #     return False

        banReact = {
            # Cu111
            # "HCOO-OadsOads.H-Hads.H2-HadsHads.H2":
            # "HCOOH-Oads.H-Hads.H-Hads.H2",

            # 2ZnCu211
            # "H2-HadsHads.CO2-CadsOads.H2.H2": "H-Hads.COOH-CadsOads.H2.H2",
            # 1ZnCu211
            # "H2COOH-Oads.H-Hads.H2": "HCHO-Oads.H2O.H2",

            # high Zn 211
            # "HCOO-OadsOads.H-Hads.H2.H2": "H2CO-OadsOads.H2.H2"

            # HCOO-1Zn211
            # "H2.H.CO2.CO2-CadsOads.H2.H2": "H2.COOH-CadsOads.CO2.H2.H2",
            # "3540832538348544.CO2.CO2-CadsOads.H2.H2":
            # "H.H2.CO2.CO2-CadsOads.H2.H2",
            # "H2.COOH-CadsOads.CO2.H2.H2":
            # "3540832538348544.CO2.CO2-CadsOads.H2.H2"

        }
        if banReact.get(name1, 0) == name2 or banReact.get(name2, 0) == name1:
            return False

        return True

    def RultOutPath(self, name1):
        for i in range(len(self.allpath) - 1, -1, -1):
            self.allpath[i].EndtoendSort(name1)
            # print(i)
            for pair in self.allpath[i][0:-1]:
                if (self.AllowedReactRule(pair[0].sminame,
                                          pair[1].sminame)) is False:
                    print("POP:" + pair[0].sminame + "--" + pair[1].sminame)
                    self.allpath.pop(i)
                    break

    def PathFinderYenAlgor(self, graph, source, slink, kpath, printFlag=False):
        # Improved from https://zhuanlan.zhihu.com/p/39709567
        A = []
        B = []
        costList = []
        path = nx.dijkstra_path(graph, source, slink)
        cost = 0
        for nodei in range(len(path) - 1):
            cost += graph.get_edge_data(path[nodei], path[nodei + 1])["weight"]
        costList.append(cost)
        A.append(tuple(path))
        for k in range(1, kpath):
            for i in range(0, len(A[k - 1]) - 1):
                g = cp.deepcopy(graph)
                spurnode = A[k - 1][i]
                rootpath = A[k - 1][0:i]
                for p in A:
                    if rootpath == p[0:i] and spurnode == p[i]:
                        if g.has_edge(p[i], p[i + 1]):
                            g.remove_edge(p[i], p[i + 1])
                for rootpathnode in rootpath:
                    if g.has_node(rootpathnode):
                        g.remove_node(rootpathnode)
                try:

                    spurpath = nx.dijkstra_path(g, spurnode, slink)
                    totalpath = list(rootpath) + spurpath

                    flag = True
                    for path in B:
                        if totalpath == path[0]:
                            flag = False
                            break
                    if flag:
                        B.append([totalpath, 0])

                except nx.NetworkXNoPath:
                    continue

            if not B:
                break
            # g = cp.deepcopy(graph)

            for path in B:
                if path[1] == 0:
                    # calculate cost
                    for nodei in range(len(path[0]) - 1):
                        path[1] += graph.get_edge_data(
                            path[0][nodei], path[0][nodei + 1])["weight"]
            B.sort(key=lambda cost: cost[1], reverse=True)
            if printFlag:
                print("%d %e" % (k, B[-1][1]))
            costList.append(B[-1][1])
            A.append(tuple(B.pop()[0]))

        # print(A)
        return A, costList


def ShortSmiName(sminame):
    adict = {}
    for ads in sminame.split("."):
        if ads in adict:
            adict[ads] += 1
        else:
            adict[ads] = 1
    alist = []
    for ads in adict:
        if adict[ads] == 1:
            alist.append("%s" % ShortFormula(ads))
        else:
            alist.append("%d*%s" % (adict[ads], ShortFormula(ads)))
    return "+".join(alist)


def ShortFormula(name):
    formulaDict = {
        "Hydrogen_Gas": "H2",
        "Carbon_Dioxide": "CO2",
        "[H]-Hads": "H*",
        "Formic_Acid": "HCOOH",
        "Water": "H2O",
        "Formaldehyde": "HCHO",
        "Hydroxy-Oads": "OH*"
    }
    if name in formulaDict:
        return formulaDict[name]
    else:
        return name


# if __name__ == "__main__":
#   test= allstr()
#   import sys
#   test.readfile(sys.argv[1])
#   test.GetAllsminame(numproc = 6,colorflag = 0)
#   test.GetAllECFPname(numproc = 6)

#  test.AddGibbsEnergy(strType="pair", filemode=2)
#  for struc in test:
#      print((struc.sminame,struc.energy- struc.trueEnergy))

# import sys
# test = AllPair()
# test.readfile(filename="test.arc", filemode=2)
# print(len(test))
# for pair in test:
#     print(pair[0].sminame, pair[0].energy)#, pair[0].zpe)
#     print(pair.TSstr.sminame, pair.TSstr.energy)#, pair.TSstr.zpe)
#     print(pair[1].sminame, pair[1].sminame)#, pair[1].zpe)
