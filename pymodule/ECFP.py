from allstr_new import allstr
from structure_new import Str, wrapCalFakebmx, wrapSegmolecular
import numpy as np
import hashlib
# import zlib
import ctypes
from ctypes import pointer
# from babelfunc import *
from rdkitInterface import glueSegStr, glueSegStr_pure, calAllName, calmass, RdkitSminame
from multiprocessing import Pool
from PeriodicTable import Eletable

readableNameDict = {
    "124853441469919219422722":
    "CO",
    "15180224853442318479815565137":
    "CO2",
    "2795242429503218":
    "H2",
    # ---------------------------------------------
    "19558913":
    "H",
    "19933032":
    "O",
    # ---------------------------------------------
    "118827371478685318633138":
    "OH",
    "1480732210601812478685319215682":
    "H2O",
    "11882737147868531863313819558913":
    "H2O",
    "13740435347868531617256838040874":
    "H3O",

    # ---------------------------------------------
    "22485344134847971478685315450431164563811859689719707065":
    "HCOO",
    "1355693159722312483157124853441318479814786853150130041556513718013593":
    "COOH",
    "1597223122374321248315712485344125895682478685315450431167677251859689719707065":
    "HCOOH",
    "2248534413484797147868531545043116456381185968971955891319707065":
    "HCOOH",
    "259722322483157146327022478685315565137163679882642965928013593":
    "C(OH)2",
    "224853442316317424786853167634021836061128442205":
    "H2CO",
    "25972231146357321786040222374322248315734786853148643421859689719707065":
    "HC(OH)2",
    "25972232248315723163174141309912443193744786853273479921817833018360611":
    "H2C(OH)2",
    "1597223124831571248534423163174146521123478685317347992176063671836061118442205":
    "H2COOH",
    # ---------------------------------------------
    "15972231142782112483157142109923478685324804582155651371695875218013593197347501974817719884318":
    "HOCOH2",
    "1232005314210992547868532480458236787907178783801957426819935132":
    "CH3OH2",
    "1597223124831571252414523163174142109921445277254786853248045821582577817347992183606111878251619592565":
    "HOCH2OH2",

    # ---------------------------------------------
    "111010631248534414786853175857571805173618322925":
    "CH=O",
    "1597223162024412483157147868531486950019422722":
    "COH",
    "13205861597223124831571279403413001095247868531758575718051736":
    "HCOH",
    "124853442478685325153926163915541761189719358688":
    "HCHO",
    "118159411232005312485344145061803478685336787907":
    "CH3O",
    "15972231604059124831571370403414152259347868532515392616391554":
    "CH2OH",
    "159722311252378123200531248315744786853160163053678790718204349":
    "CH3OH",
    "11815941123200531248534414506180347868533678790719558913":
    "CH3OH",
    # ---------------------------------------------
    '13184126':
    "C",
    "139511971440439014786853":
    "CH",
    "24786853259944191622691418486911":
    "CH2",
    "12834414337978743478685316389191":
    "CH3",
    "4995519113345124478685318879641":
    "CH4",

    # -------------
    "1302479114210992447868532480458225153926159288311639155419603628":
    "H2COH2",
    "13225686240339861540832527952424":
    "H3",
    "180735911316635142109923478685324804582150228541758575718051736":
    "HCOH2",

    # ----------C2 species
}


def wrapGenECFPname(datain):
    str, savefp = datain
    # print str.sminamef
    ECFPname, FP = str.GenECFPname()
    if savefp:
        return ECFPname, FP
    else:
        return ECFPname


def wrapGenECFPname_Conly(datain):
    str, savefp = datain
    #print str.sminame
    ECFPname, FP = str.GenECFPname(LConly=1)
    if savefp:
        return ECFPname, FP
    else:
        return ECFPname


def one_of_k_encoding(x, allowable_set):
    if x not in allowable_set:
        raise Exception("input {0} not in allowable set{1}:".format(
            x, allowable_set))
    return list(map(lambda s: x == s, allowable_set))


def one_of_k_encoding_unk(x, allowable_set):
    """Maps inputs not in the allowable set to the last element."""
    if x not in allowable_set:
        x = allowable_set[-1]
    return list(map(lambda s: x == s, allowable_set))


def GetAllSminame_FromEcfp_Single(struc):
    fragments = [[] for i in np.unique(struc.group)]
    for id, atom in enumerate(struc.atom):
        if (atom.ele > 18):
            continue
        atom.id = id
        fragments[struc.group[id] - 1].append(id)

    _tmp = AllStr()
    for frag in fragments:
        if (len(frag) == 0):
            continue
        _tmp.append(struc.GenSubStr(frag))

    for _struc in _tmp:
        surfaceflag = "-"
        for atom in _struc.atom:
            if struc.surfaceatom[atom.id] == 1:
                surfaceflag += '%sads' % (Eletable[atom.ele - 1])
        _struc.surfaceflag = surfaceflag
    # _tmp.GetTmpFakebmx(numproc=numproc)
    # _tmp.calAllFakeBondMaxtrix(numproc=numproc)
    for _struc in _tmp:
        r = wrapCalFakebmx(_struc)
        _struc.bmx2D = np.array(r[1]).reshape(_struc.natom, _struc.natom)
        _struc.bmx1D = r[1]
        _struc.bondneed = r[2]
        _struc.Lminstr = r[0]
        _struc.surfaceatom = r[3]

    # _tmp.calAllSegmolecular(numproc=numproc)
    for _struc in _tmp:
        _struc.group = wrapSegmolecular(_struc)

    # _tmp.GetAllECFPname(numproc=numproc)
    for _struc in _tmp:
        _struc.ECFPname = wrapGenECFPname((_struc, False))

    for _struc in _tmp:
        _struc.name = readableNameDict.get(_struc.ECFPname, _struc.ECFPname)

        # if name could not be found
        if _struc.name is _struc.ECFPname:
            # _struc.name = RdkitSminame(_struc.atom, _struc.bmx2D, _struc.bondneed, 2,
            #              _struc.surfaceatom)
            _nameDictArc = allstr()
            _nameDictArc.append(_struc)
            _nameDictArc.printall(_struc.ECFPname + ".arc")

        if len(_struc.surfaceflag) != 1:
            _struc.name += _struc.surfaceflag
    struc.sminame = ".".join([struc.name for struc in _tmp])

    return struc.sminame


def GetAllSminame_FromEcfp_SingleProc(
    struc,
    flag=2,
):
    if flag == 1:
        bmx = struc.Bondmatrix()
        struc.bmx2D = np.array(bmx).reshape(struc.natom, struc.natom)
        struc.bmx1D = bmx
    if flag == 2:
        r = struc.CheckMin()
        struc.bmx2D = np.array(r[1]).reshape(struc.natom, struc.natom)
        struc.bmx1D = r[1]
        struc.bondneed = r[2]
        struc.Lminstr = r[0]
        struc.surfaceatom = r[3]
    struc.group = struc.Segmolecular()
    GetAllSminame_FromEcfp_Single(struc)

class SurfaceStr(Str):
    def GenFeature(self):
        self.AddAtomID()
        self.GetAtomInfo()
        # print(self.surfaceatom)
        for iatom, atom in enumerate(self.atom):
            atom.surface = self.surfaceatom[iatom]
            atom.surfacemetal = self.atom[-1].elesymbol
            atom.feature = self.Atom_Features(atom)
        self.NeighbourList()

    def Atom_Features(self, atom):
        # return np.array(one_of_k_encoding_unk(atom.elesymbol,['C', 'O', 'H', 'Cu', 'Unknown']) +
        return np.array(
            one_of_k_encoding_unk(
                atom.elesymbol,
                ['C', 'O', 'H', 'N', 'Cu', 'Pt', 'Rh', 'Zn', 'Unknown']) +
            one_of_k_encoding_unk(
                atom.expbond,
                [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 'Unknown']) +
            one_of_k_encoding_unk(atom.imph, [0, 1, 2, 3, 4, 5, 'Unknown']) +
            one_of_k_encoding_unk(atom.surfacemetal,
                                  ['Cu', 'Pt', 'Rh', 'Zn', 'Unknown']) +
            # one_of_k_encoding_unk(atom.value, [-2, -1, 0, 1, 2, 3])+
            # one_of_k_encoding_unk(atom.surface, ['single', 'bridge','hcp', 'fcc','4','unkonwn'])+
            [0])

    def Bond_Features(self, bondorder):
        return np.array(one_of_k_encoding(bondorder, [1, 2, 3, 4, 5]))


#    def GetNeighbour(self,iatom):
#        for i in range(self.natom):
#            if self.bmx2D[iatom][i]>0:
#                self.atom[iatom].neighbour.append()
#

    def NeighbourList(self):
        for i in range(self.natom):
            self.atom[i].nbatom = []
            self.atom[i].bond = []
            self.atom[i].nb = []

        for i in range(self.natom):
            for j in range(i + 1, self.natom):
                if self.bmx2D[i][j] > 0:
                    # do i need skip H?
                    self.atom[j].nbatom.append(self.atom[i])
                    self.atom[i].nbatom.append(self.atom[j])
                    # self.atom[j].bond.append(self.Bond_Features(self.bmx2D[i][j]))
                    # self.atom[i].bond.append(self.Bond_Features(self.bmx2D[i][j]))
                    self.atom[j].nb.append([self.bmx2D[i][j], i])
                    self.atom[i].nb.append([self.bmx2D[i][j], j])

    def SortbyDegree(self, ):
        self.atom.sort(key=lambda x: len(x.neighbour))

    def AllAtomFeature(self):
        atomfeature = []
        for atom in self.atom:
            atomfeature.append(atom.feature)
            atomfeature.append(atom.bond)

    # def AllBondFeature(self):
    #     atomfeature = []
    #     for atom in self.atom:
    #         for bond in atom.bond:
    #             bondfeature.append(bond)
    #             bondfeature.append(bond)

    def GenECFPname(self, depth=3, dim=10000000, savefp=False, LConly=0):
        FP = ECFP(self, depth, dim)
        FP.GenECFP()
        string = ''
        for index in FP.allindex:
            if (LConly and FP.IndextoFp[index][0].coreele != 6):
                continue
            # print index
            string = string + str(len(FP.IndextoFp[index]))
            string = string + str(index)
        self.ECFPname = string
        if savefp:
            self.FP = FP

        # print self.ECFPname
        return string, FP

    def GenSubStr(self, atomlist):
        substr = SurfaceStr()
        substr.atom = []
        for i in atomlist:
            substr.atom.append(self.atom[i])
        substr.abc = self.abc
        substr.energy = 0.0
        substr.sortatombyele()
        substr.calAtomnum()
        substr.abc2lat()
        return substr

    def CheckMin(self, flag=1):
        sqnatm = self.natom**2

        self.calCtypes()
        # print("use Lib_fillbond_tmp")
        # program='/home10/kpl/pymodule/Lib/Lib_fillbond_tmp/checkminbond.so'
        program = "/home9/shiyf/bin/kpl-path-shiyfgit/Lib_fillbond_tmp/checkminbond.so"
        # program='/home10/kpl/pymodule/Lib/Lib_fillbond/checkminbond.so'
        Lminstr = pointer(ctypes.c_bool(0))
        bmatrix = pointer((ctypes.c_int * sqnatm)(*[0 for i in range(sqnatm)]))
        bondneed = pointer(
            (ctypes.c_int * (self.natom))(*[0 for i in range(self.natom)]))
        surface = pointer(
            (ctypes.c_int * (self.natom))(*[0 for i in range(self.natom)]))

        checkmin = ctypes.cdll.LoadLibrary(program)
        if flag == 0:
            checkmin.judgebond_(self.c_natm, self.c_iza, self.c_xa, self.c_rv,
                                Lminstr, bmatrix, bondneed)
        elif flag == 1:
            # print 'into surface'
            checkmin.judgebondsurface_(self.c_natm, self.c_iza, self.c_xa,
                                       self.c_rv, Lminstr, bmatrix, bondneed,
                                       surface)
            # print 'finish fffffffffffffff'
        # print(flag)
        bmx = list(bmatrix.contents)
        # print(program)
        # self.bmx2D = np.array(bmx).reshape(self.natom, self.natom)
        # self.bmx1D = bmx
        #    bmx2D=np.array(bmx).reshape(self.natom, self.natom)
        #    for i in range(0,9):
        #         for j in range(i+1,9):
        #            #print(i,j)
        #            if abs(bmx2D[i,j]) > 0.01 :
        #                print(i+1,j+1)
        #    print("1")
        self.bondneed = list(bondneed.contents)
        self.Lminstr = bool(Lminstr.contents)
        return self.Lminstr, bmx, list(bondneed.contents), list(
            surface.contents)


class AllStr(allstr):
    def readfile(self, inputfile, forcefile=False, allformat=0):
        f = open(inputfile, 'r')
        currentStr = -1
        for line in f:
            if ('Energy' in line or 'React' in line
                    or ('TS' in line and 'sect' not in line)
                    or 'SSW-Crystl' in line or 'Str' in line):
                self.append(SurfaceStr())
                currentStr += 1
                self[currentStr].Lfor = False
                try:
                    self[currentStr].energy = float(line.split()[-1])
                    try:
                        self[currentStr].maxFF = float(line.split()[-2])
                    except Exception:
                        self[currentStr].maxFF = 0

                    if self[currentStr].energy.is_integer():
                        self[currentStr].energy = float(line.split()[-2])
                except Exception:
                    self[currentStr].energy = float(line.split()[-2])
                    self[currentStr].maxFF = 0
            elif 'CORE' in line:
                self[currentStr].addatom(line, 1)
            elif ('PBC' in line) and ('ON' not in line):
                self[currentStr].abc = [float(x) for x in line.split()[1:]]
        f.close()
        for str in self:
            str.sortatombyele()
            str.calAtomnum()
            str.abc2lat()

        if forcefile:
            f = open(forcefile, 'r')
            currentStr = -1
            for line in f:
                if "For" in line:
                    self[currentStr].Lfor = True
                    currentStr += 1
                    iatom = 0
                    for atom in self[currentStr].atom:
                        atom.force = [0.0, 0.0, 0.0]
                elif len(line.split()) == 6:
                    self[currentStr].addStress(line)
                elif len(line.split()) == 3:
                    if "****" not in line:
                        self[currentStr].addForce(line, iatom)
                    else:
                        self[currentStr].addForce('0.0 0.0 0.0', iatom)
                    iatom += 1

        if allformat:
            for str in self:
                str.TransferToXYZcoordStr()

    def GetAllECFPname(self, numproc=24, savefp=False, LConly=False):
        _tmp = []
        for struc in self:
            _tmp.append((struc, savefp))

        pool = Pool(processes=numproc)
        if LConly:
            result = pool.map_async(wrapGenECFPname_Conly, _tmp)
        else:
            result = pool.map_async(wrapGenECFPname, _tmp)
        pool.close()
        pool.join()

        for struc, r in zip(self, result.get()):
            if savefp:
                struc.ECFPname = r[0]
                struc.FP = r[1]
            else:
                struc.ECFPname = r

    def GetAllFrag(self, numproc=24, flag=2):
        if flag == 1:
            self.calAllBondMatrix(numproc=numproc)
        if flag == 2:
            self.calAllFakeBondMaxtrix(numproc=numproc)
        self.calAllSegmolecular(numproc=numproc)

        outstr = AllStr()
        for struc in self:
            fragments = [[] for i in np.unique(struc.group)]
            for id, atom in enumerate(struc.atom):
                if (atom.ele > 18):
                    continue
                atom.id = id
                fragments[struc.group[id] - 1].append(id)
            for frag in fragments:
                if (len(frag) == 0):
                    continue
                outstr.append(struc.GenSubStr(frag))

        return outstr

    def GetTmpFakebmx(self, numproc=24, flag=2, colorflag=1):
        if flag == 1:
            self.calAllBondMatrix(numproc=numproc)
        if flag == 2:
            self.calAllFakeBondMaxtrix(numproc=numproc)
        self.calAllSegmolecular(numproc=numproc)

    def GetAllSminame_FromEcfp(self,
                               numproc=24,
                               flag=2,
                               colorflag=1,
                               multiProc=True):
        if flag == 1:
            self.calAllBondMatrix(numproc=numproc)
        if flag == 2:
            # print("tesst")
            self.calAllFakeBondMaxtrix(numproc=numproc)
        self.calAllSegmolecular(numproc=numproc)

        if multiProc is True:
            pool = Pool(processes=numproc)
            result = pool.map_async(GetAllSminame_FromEcfp_Single, self)
            pool.close()
            pool.join()
            for struc, sminame in zip(self, result.get()):
                struc.sminame = sminame
        else:
            for struc in self:
                struc.sminame = GetAllSminame_FromEcfp_Single(struc)

        # for struc in self:
        #     GetAllSminame_FromEcfp_Single(struc)

    def GetAllsminame_fromreactstr(self, numproc=24, flag=2, colorflag=1):
        for struc in self:
            struc.RemoveMetalBond()
        self.calAllSegmolecular(numproc=numproc)

        allgroup = []
        for str in self:
            substr = [[] for i in np.unique(str.group)]
            for id, atom in enumerate(str.atom):
                atom.id = id
                substr[str.group[id] - 1].append(atom)
            substr = sorted(substr, key=lambda x: calmass(x), reverse=True)
            if flag == 1:
                allgroup.append((substr, str.bmx2D, str.lat, [], 1, []))
            if flag == 2:
                allgroup.append((substr, str.bmx2D, str.lat, str.bondneed, 2,
                                 str.surfaceatom))
            # print str.bondneed

        pool = Pool(processes=numproc)
        result = pool.map_async(calAllName, allgroup)
        pool.close()
        pool.join()

        for istr, (str, re) in enumerate(zip(self, result.get())):
            str.allmol = re
            str.id = istr

        for str in self:
            if colorflag:
                str.sminame, strflag = glueSegStr(str.allmol)
            else:
                str.sminame, strflag = glueSegStr_pure(str.allmol)

        for str in self:
            str.bmx2D = str.bmxsave


class singlefp(object):
    def __init__(self):
        self.index = 0
        self.neighbour = []
        self.allatom = []
        self.trace = []
        self.save = 1
        self.coreele = 0
        self.coreid = 0

    def Inherit(self, parent):
        self.index = parent.index
        self.allatom = parent.allatom[:]
        self.trace = parent.trace[:]
        self.coreele = parent.coreele
        self.coreid = parent.coreid

    def genidentifier(self, dim=10000000):
        self.neighbour.sort(key=lambda x: x[0])

        # if self.coreele == 29:
        #    print len(self.neighbour)

        for nb in self.neighbour:
            self.trace.append(nb[1])

        self.allatom.sort()
        self.allatomset = np.unique(self.allatom)

        tmpall = [(1, self.index)]
        for nb in self.neighbour:
            tmpall.append(nb[0])
            # print nb[0]
            # tmpall.append(nb[0][0])
            # tmpall.append(nb[0][1])
        self.index = hashlist(tmpall, dim)

        # string = '1:%s'(%self.index)
        # for item in self.neighbour:
        #    string = string+ str(item)


class ECFP(object):
    def __init__(self, structure, nlayer, dim=10000000):
        self.str = structure
        self.nlayer = nlayer
        self.allfp = []
        self.dim = dim

#    def PreparebyLayer(self):
#        for atom in self.str.atom:
#            atom.layeratom = {}
#            atom.layeratom[0] = [atom]
#
#        neighbourdict = {}
#        for ilayer in range(nlayer):
#            atom.layeratom[ilayer+1] = []
#            for iatom,atom in enumerate(self.str.atom):
#                #neighbourdict[(iatom,ilayer)] = []
#                for subatom in atom.layeratom[ilayer]:
#                    #neighbourdict[(iatom,ilayer)].append(subatom.id)
#                    for inb,newatom in enumerate(subatom.neighbour):
#                        atom.layeratom[ilayer+1].append(newatom)
#                        self.cyclenb[(iatom,ilayer+1)].append([subatom.bond[inb],newatom.id])

    def GenECFP(self):
        self.InitIdentifier()
        for ilayer in range(1, self.nlayer):
            self.UpdateLayer(ilayer)
        allindex = []
        for fp in self.allfp:
            allindex.append(fp.index)
        allindex.sort()
        self.allindex = np.unique(allindex)
        # self.allindex,self.allindexconut =list(np.unique(allindex,return_counts=True))
        self.GenDict()
        return

    def GenDict(self):
        # AtomsettoIndex = {}
        IndextoAtomset = {}
        IndextoFp = {}
        for fp in self.allfp:
            # AtomsettoIndex[set(fp.allatomset)] = fp.index
            IndextoAtomset[fp.index] = fp.allatomset
            if fp.index not in IndextoFp.keys():
                IndextoFp[fp.index] = [fp]
            else:
                IndextoFp[fp.index].append(fp)

        self.IndextoAtomset = IndextoAtomset
        self.IndextoFp = IndextoFp

    def InitIdentifier(self):
        self.str.GenFeature()
        self.genbylayer = [[0 for x in range(self.str.natom)]]
        for iatom, atom in enumerate(self.str.atom):
            newfp = singlefp()
            newfp.allatom.append(atom.id)
            newfp.trace.append(atom.id)
            newfp.coreele = atom.ele
            newfp.coreid = iatom
            # print atom.feature
            newfp.index = hashlist(atom.feature, self.dim)
            # print newfp.index
            # print newfp.index
            newfp.allatom.sort()
            newfp.allatomset = np.unique(newfp.allatom)
            # newfplist.append(newfp)
            self.genbylayer[0][iatom] = newfp
            # if newfp.coreele < 18: continue
            # if newfp.coreele == 1: continue
            self.allfp.append(newfp)

    def UpdateLayer(self, layer):
        newfplist = []
        for iatom, atom in enumerate(self.str.atom):
            newfp = singlefp()
            # newfp.append((1,self.genbylayer[layer-1][atom.id]))
            newfp.Inherit(self.genbylayer[layer - 1][atom.id])
            for subatomid in self.genbylayer[layer - 1][atom.id].allatom:
                # if self.str.atom[subatomid].ele > 18: continue
                for nb in self.str.atom[subatomid].nb:
                    if nb[1] not in newfp.allatom:  # and self.str.atom[nb[1]].ele != 1:
                        newfp.neighbour.append([
                            (nb[0], self.genbylayer[layer - 1][nb[1]].index),
                            nb[1]
                        ])
                        newfp.allatom.append(nb[1])

            # newfp.allatomset = set(newfp.allatom)
            newfp.genidentifier(self.dim)
            newfplist.append(newfp)
        # for fp in  newfplist:
        # print fp.save, fp.index
        newfplist = self.CheckDuplicate(newfplist)
        # print 'CheckDuplicate'
        # for fp in  newfplist:
        #    print fp.save, fp.index

        self.genbylayer.append([0 for x in range(self.str.natom)])
        for ifp, fp in enumerate(newfplist):
            if fp.save == 1:
                self.genbylayer[layer][ifp] = fp
                # if fp.coreele < 18: continue
                # if fp.coreele == 1: continue
                self.allfp.append(fp)
            else:
                self.genbylayer[layer][ifp] = self.genbylayer[layer - 1][ifp]
        return

    def CheckDuplicate(self, newfplist):
        for newfp in newfplist:
            for fp in self.allfp:
                if set(newfp.allatomset) == set(fp.allatomset):
                    # for unity
                    # newfp.index = fp.index
                    # old mode
                    newfp.save = 0
                    break

        for newfp1 in newfplist:
            for newfp2 in newfplist:
                if set(newfp1.allatomset) == set(newfp2.allatomset):
                    if newfp1.index > newfp2.index:
                        # newfp1.index = newfp2.index
                        newfp1.save = 0
                    elif newfp1.index < newfp2.index:
                        # newfp2.index = newfp1.index
                        newfp2.save = 0
        return newfplist


# class NNFingerPrint(object):
#    def __init__(self):
#        return
#
#    def neighbouractive(self):
#        for degree in degrees:
#            atom_neighbour= self.atom[j].neighbour
#            bond_neighbour = self.bond[j].neighbour
#
#
#    def layerupdate(self):
#        return


def hashlist(list, dim=10000000):
    # print list
    string = ''
    for i in list:
        string = string + str(i)

    # print string
    md5 = hashlib.md5()

    # print md5.digest_size
    # print md5.block_size
    md5.update(string.encode())
    # print md5.digest
    # print string

    # return md5.hexdigest()
    _tmp = md5.hexdigest()
    # print int((int(_tmp,16))%dim)
    return int((int(_tmp, 16)) % dim)
    # return zlib.adler32(string)


if (__name__ == "__main__"):
    test = [1, 2332421, 1, 14757242, 2, 124253533]
    #    string = ''
    #    for i in test:
    #        string = string +str(i)
    #    sha1 =hashlib.sha1()
    #    sha1.update(string)
    #    print sha1.hexdigest()
    #
    #
    #
    #    md5 =hashlib.md5()
    #    md5.update(string)
    #    print md5.hexdigest()
    #
    #
    #
    #    print zlib.adler32(string)
    print(hashlist(test))

    test = AllStr()
    test.readfile('target.arc')
    test.calAllFakeBondMaxtrix()
    fp = ECFP(test[0], 3)
    fp.GenECFP()
    # print len(fp.allfp)

    # for sfp in fp.allfp:
    #    print sfp.index

#    print (fp.allindex)
#    print (fp.IndextoAtomset)
# for ilayer in range(3):
#     for iatom in range(test[0].natom):
#         print (ilayer,iatom,fp.genbylayer[ilayer][iatom].index, fp.genbylayer[ilayer][iatom].trace)
#    all.GetAllpuresminame()
#    all.GenCanoicalSmiles(chiral= False)
#
#    ndim = 1000000
#    info= {}
#    info2 = {}
#    for i in range(len(test)):
#        str1 = test[i]
#        m1 = Chem.MolFromSmiles(str1.canname)
#        Chem.SanitizeMol(m1)
#        fptest= AllChem.GetHashedMorganFingerprint(m1,2,nBits=ndim, bitInfo = info,useChirality=False)
#        if i ==0:
#            fpall =fptest
#        fpall =fpall + fptest
#    print info
