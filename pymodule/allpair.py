# from allstr_uncm import AllStr
from ECFP import AllStr
# import os
# import numpy as np
import ctypes
from ctypes import pointer
from multiprocessing import Pool


def wrapCalISFSdistance(datain):
    pair = datain
    # print str.sminame
    dis = pair.ISFSdistance()
    return dis


class ReactionPair(AllStr):
    def __init__(self, str1, str2, TS=0, TSstr=False, sortflag=True):
        list.__init__(self)
        self.append(str1)
        self.append(str2)
        # self.append(TSstr)
        self.TS = TS
        self.TSstr = TSstr

        if sortflag:
            self.sort(key=lambda X: X.name)
        self.name = '%s' % (self[0].name + self[1].name)
        return

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


class AllPair(AllStr):
    def readfile(self, filename, Lscreenupper=0):
        _tmp = AllStr()
        _tmp.readfile(filename)
        if Lscreenupper == 1:
            for Str in _tmp:
                Str.screenuppersurf()

        _tmp.GetAllSminame_FromEcfp(numproc=8, colorflag=0)
        _tmp.GetAllECFPname()

        for i in range(0, len(_tmp), 2):
            if (not _tmp[i].sminame) or (not _tmp[i + 1].sminame):
                continue
            if Lscreenupper == 1:
                if _tmp[i].upper == 1 or _tmp[i + 1].upper == 1:
                    continue
            _tmp[i].name = _tmp[i].ECFPname
            _tmp[i + 1].name = _tmp[i + 1].ECFPname
            self.append(ReactionPair(_tmp[i], _tmp[i + 1]))


#    def GetAlldistance(self):
#        for pair in self:
#            pair.ISFSdistance()

    def GetAlldistance(self, numproc=24):
        _tmp = []
        for pair in self:
            _tmp.append(pair)

        pool = Pool(processes=numproc)
        result = pool.map_async(wrapCalISFSdistance, _tmp)
        pool.close()
        pool.join()

        for pair, r in zip(self, result.get()):
            pair.totdis = r

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
        self.allstrmin = {}
        for name in self.allname.keys():
            self.allstrmin[name] = min(self.allstrbyname[name],
                                       key=lambda x: x.energy)

    def GenPairDict(self):
        # self.GetAlldistance()
        self.allpair = {}
        for pair in self:
            if pair.name not in self.allpair.keys():
                self.allpair[pair.name] = []
                self.allpair[pair.name].append(pair)
            else:
                self.allpair[pair.name].append(pair)

        # self.shortpair = {}
        # for name in self.allpair.keys():
        #     self.allpair[name].sort(key=lambda x: x.totdis)
        #     self.shortpair[name] = min(self.allpair[name],
        #                                key=lambda x: x.totdis)

    def OutPair(self, nmaxtotal=1000, nmax=50, outfile='outpair.arc'):
        _tmp = AllStr()
        if len(self) < nmaxtotal:
            for pair in self:
                _tmp.append(pair[0])
                _tmp.append(pair[1])
        else:
            self.GenPairDict()
            for name in self.allpair.keys():
                # print name, len(self.allpair[name])
                if len(self.allpair[name]) > nmax:
                    for pair in self.allpair[name][:nmax]:
                        _tmp.append(pair[0])
                        _tmp.append(pair[1])
                else:
                    for pair in self.allpair[name]:
                        _tmp.append(pair[0])
                        _tmp.append(pair[1])
        _tmp.printall(outfile)
        return

    def GenPairdict_byname(self):
        self.AllName()
        self.GenPairDict()
        self.rpdict = {}
        for name in self.allname.keys():
            self.rpdict[name] = []
            for pair in self.shortpair.values():
                if pair[0].name == name or pair[1].name == name:
                    self.rpdict[name].append(pair)
            self.rpdict[name].sort(key=lambda X: X.totdis)
