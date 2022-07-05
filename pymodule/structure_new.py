from atom_k import S_atom
from bond_k import bond
import numpy as np
import PeriodicTable as PT
import math as m
import os
import ctypes
from ctypes import pointer
import re
from functools import reduce


def wrapCalBondMatrix(str):
    return str.Bondmatrix()


def wrapSegmolecular(str):
    return str.Segmolecular()


def wrapCalFakebmx(str):
    Lminstr, bmx, bondneed, surface = str.CheckMin()
    return Lminstr, bmx, bondneed, surface


def wrapCalQ(str):
    return str.cal_Q()


class Str(object):
    def __init__(self):
        self.atom = []
        self.lat = []
        self.natom = 0
        self.nele = 0
        self.energy = 0
        self.serial_num = 0
        self.eleList = []
        self.natompe = []
        self.stress = []
        self.abc = []
        self.Lfor = True
        self.MaxFF = 0

        self.Nat = 0
        self.EleNam = []
        self.iza = []
        self.Coord = []
        self.sp = {}
        self.sporder = {}
        self.frac = []
        self.centerf = {}
        self._cycle = []
        self.cart = []
        self.For = []

    def addatom(self, line, flag=1):
        if flag == 1:
            # for arc and ms
            coord = [float(x) for x in line.split()[1:4]]
            elesymbol = line.split()[0]
            ele = PT.Eledict[elesymbol]
            satom = S_atom(coord, ele)
            try:
                satom.charge = float(line.split()[-2])
            except Exception:
                satom.charge = 0.0
            self.atom.append(satom)
            return

        elif flag == 2:
            coord = [float(x) for x in line.split()[2:5]]
            ele = int(line.split()[1])
            self.atom.append(S_atom(coord, ele))
        elif flag == 3:  # for cat file
            coord = [float(x) for x in line.split()[1:4]]
            elesymbol = str(line.split()[0])[0]
            ele = PT.Eledict[elesymbol]
            self.atom.append(S_atom(coord, ele))
        elif flag == 4:  # for mol file
            coord = [float(x) for x in line.split()[0:3]]
            elesymbol = line.split()[3]
            ele = PT.Eledict[elesymbol]
            self.atom.append(S_atom(coord, ele))
        elif flag == 5:  # for QM9
            coord = [float(x) for x in line.split()[1:4]]
            elesymbol = line.split()[0]
            ele = PT.Eledict[elesymbol]
            self.atom.append(S_atom(coord, ele))

        elif flag == 6:
            coord = [float(x) for x in line.split()[1:4]]
            elesymbol = line.split()[-2]
            ele = PT.Eledict[elesymbol]
            satom = S_atom(coord, ele)
            satom.charge = float(line.split()[-1])
            self.atom.append(satom)

    def addForce(self, line, serial_num, flag=1):
        if flag == 1:
            # for arc
            self.atom[serial_num].force = [float(x) for x in line.split()]
        elif flag == 2:
            self.atom[serial_num].force = [float(x) for x in line.split()[2:5]]

    # def addCharge(self, base):
    #     for atom in self.atom:
    #         atom.charge = base[atom.elesymbol]
    def addCharge(self):
        for atom in self.atom:
            atom.charge = 0

    def GetMaxF(self):
        maxf = 0
        for atom in self.atom:
            _tmpf = max([abs(x) for x in atom.force])
            if maxf < _tmpf:
                maxf = _tmpf
        self.maxf = maxf

    def addStress(self,
                  line,
                  flag=1):  # type 1 => for arc file; 2 => Data file
        if flag == 1:
            self.stress = [float(x) for x in line.split()]
        elif flag == 2:
            self.stress = [float(x) for x in line.split()[1:7]]

    def calAtomnum(self, ):
        self.natom = len(self.atom)
        self.eleList, self.natompe = list(
            np.unique([atom.ele for atom in self.atom], return_counts=True))
        self.nele = len(self.eleList)
        self.elenameList = [PT.Eletable[ele - 1] for ele in self.eleList]
        self.eleCompos = dict(zip(self.eleList, self.natompe))
        self.eleDict = self.eleCompos

    def screenuppersurf(self, ):
        self.upper = 0
        for atom in self.atom:
            if atom.xyz[2] > 0.9 * self.abc[2]:
                self.upper = 1
                break

    def sortatombyele(self, ):
        self.atom.sort(key=lambda X: X.ele)

    def sortatombyz(self, ):
        self.atom.sort(key=lambda X: X.xyz[2])

    def calTwoDimCoord(self, ):
        self.xa = np.array([atom.xyz for atom in self.atom])

    def calOneDimCoord(self, ):
        self.calTwoDimCoord()
        self.xa_onedim = reduce(lambda a, b: a + b, [list(x) for x in self.xa])

    def AddAtomID(self):
        for iatom, atom in enumerate(self.atom):
            atom.id = iatom

    def cdnt2fcnt(self, ):
        self.calTwoDimCoord()
        latinv = np.linalg.inv(self.lat)
        self.fdnt = [list(x) for x in np.matmul(self.xa, latinv)]

    def fdnt2cdnt(self, ):
        fa = np.array(self.fdnt)
        self.xa = np.matmul(fa, self.lat)
        for at, cdnt in zip(self.atom, self.xa):
            at.xyz = list(cdnt)

    def cal_dis(self, iatom1, iatom2):
        self.cdnt2fcnt()
        vbond = np.array(self.fdnt[iatom1]) - np.array(self.fdnt[iatom2])
        dis = np.linalg.norm(
            np.matmul(np.array([x - np.round(x) for x in vbond]), self.lat))
        return dis

    def cal_dis_matrix(self):
        self.Dmatrix = np.zeros([self.natom, self.natom])
        for i in range(self.natom):
            for j in range(i + 1, self.natom):
                self.Dmatrix[i, j] = self.cal_dis(i, j)
                self.Dmatrix[j, i] = self.Dmatrix[i, j]
        return

    def calcentroid(self):
        self.calAtomnum()
        xyzall = np.array([0, 0, 0])
        for atom in self.atom:
            xyzall = xyzall + np.array(atom.xyz)
        self.centroid = xyzall / self.natom

    def mvstrtoboxcenter(self):
        # for Ortholat
        self.calTwoDimCoord()
        self.calcentroid()
        center = np.array([self.abc[0] / 2, self.abc[1] / 2, self.abc[2] / 2])
        mv = center - self.centroid
        for i, atom in enumerate(self.atom):
            atom.xyz = list(self.xa[i] + mv)

    def addOrtholat(self):
        # self.calAtomnum()
        if self.natom != 0:
            self.calTwoDimCoord()
            max3d = np.amax(self.xa, axis=0)
            min3d = np.amin(self.xa, axis=0)
            delta = max(list(max3d - min3d))
            a = 10
            while ((a - delta) < 4):
                a = a + 5
            self.abc = [a, a, a, 90, 90, 90]

    def abc2lat(self):
        a, b, c = self.abc[0:3]
        alpha, beta, gamma = [x * np.pi / 180.0 for x in self.abc[3:]]

        bc2 = b**2 + c**2 - 2 * b * c * m.cos(alpha)
        h1 = a
        # h1= checkzero(h1)
        h2 = b * m.cos(gamma)
        # h2 = checkzero(h2)
        h3 = b * m.sin(gamma)
        #  h3 = checkzero(h3)
        h4 = c * m.cos(beta)
        # h4 = checkzero(h4)
        h5 = ((h2 - h4)**2 + h3**2 + c**2 - h4**2 - bc2) / (2 * h3)
        # h5 = checkzero(h5)
        h6 = m.sqrt(c**2 - h4**2 - h5**2)
        # h6 = checkzero(h6)

        self.lat = [[h1, 0., 0.], [h2, h3, 0.], [h4, h5, h6]]

    def lat2abc(self, flag=False, inlat=False):
        if not flag:
            lat = self.lat
        else:
            lat = inlat

        self.nlat = np.array(self.lat)
        a = np.linalg.norm(lat[0])
        b = np.linalg.norm(lat[1])
        c = np.linalg.norm(lat[2])
        alpha = m.acos(np.dot(lat[1], lat[2]) / (b * c)) * 180.0 / np.pi
        beta = m.acos(np.dot(lat[0], lat[2]) / (a * c)) * 180.0 / np.pi
        gamma = m.acos(np.dot(lat[0], lat[1]) / (a * b)) * 180.0 / np.pi
        if not flag:
            self.abc = [a, b, c, alpha, beta, gamma]
        else:
            return [a, b, c, alpha, beta, gamma]

    def outPOSCAR(self, outfile):
        f = open(outfile, 'w')
        f.write('system\n')
        f.write('1.000000000000\n')
        for item in self.lat:
            f.write('%12.8f  %12.8f  %12.8f\n' % (item[0], item[1], item[2]))
        f.write(
            "%s\n" %
            reduce(lambda a, b: a + b, ["%4s" % s for s in self.elenameList]))
        f.write(
            "%s\n" %
            reduce(lambda a, b: a + b, ["%4d" % num for num in self.natompe]))
        f.write('Cart\n')
        for atom in self.atom:
            f.write('%12.8f  %12.8f  %12.8f\n' %
                    (atom.xyz[0], atom.xyz[1], atom.xyz[2]))
        f.close()

    def genPOTCAR(self, sourcedir, outfile):
        os.system('rm -rf %s' % outfile)
        for ele in self.elenameList:
            os.system('cat %s/POTCAR.%s >> %s' % (sourcedir, ele, outfile))

#    def cal_distance(self,iatom1,iatom2):
#        d2 =0
#        for i in range(3):
#            d2 =d2 +(self.atom[iatom1].xyz[i]-self.atom[iatom2].xyz[i])**2
#        dis= m.sqrt(d2)
#        return dis

    def cal_neighbour(self, iatom):
        dict = {}
        for i in range(self.natom):
            if (i != iatom):
                dis = self.cal_dis(iatom, i)
                dict[i] = dis
        return dict

    def specialneighbour(self, iatom, spe_ele):
        dict = {}
        for i in range(self.natom):
            if (i != iatom and self.atom[i].ele == spe_ele):
                dis = self.cal_dis(iatom, i)
                dict[i] = dis
        return dict

    def longestbond(self, ele1, ele2, lim):
        self.cdnt2fcnt()
        maxlist = []
        for i, atom in enumerate(self.atom):
            if atom.ele == ele1:
                result = self.specialneighbour(i, ele2)
                dis = [x for x in result.values() if x < lim]
                if len(dis) > 0:
                    maxlist.append(max(dis))
        longest = max(maxlist)
        return longest

    def shortestbond(self, ele1, ele2, lim):
        self.cdnt2fcnt()
        minlist = []
        for i, atom in enumerate(self.atom):
            if atom.ele == ele1:
                result = self.specialneighbour(i, ele2)
                dis = [x for x in result.values()]
                if len(dis) > 0:
                    minlist.append(min(dis))
        short = min(minlist)
        return short

    def simpleclass(self, iatom):
        for i in range(self.natom):
            if (i != iatom):
                # dislist.append(self.cal_distance(iatom,i))
                # print dis
                dis2 = self.cal_dis(iatom, i)
                # print dis
                bondtest = bond(self.atom[iatom].ele, self.atom[i].ele, dis2)
                bondorder = bondtest.judge_bondorder()
                if bondorder != 0:
                    if self.atom[i].ele == 8:
                        self.atom[iatom].bondtype.append(101)
                    if self.atom[i].ele == 6:
                        self.atom[iatom].bondtype.append(11)
                    if self.atom[i].ele == 1:
                        self.atom[iatom].bondtype.append(1)

        self.atom[iatom].bondtype.sort()
        return

    def bondsearch(self, iatom, ele2, flag=1):
        bondlist = []
        for i in range(self.natom):
            if (i != iatom) and self.atom[i].ele == ele2:
                if flag == 2 and i < iatom and self.atom[iatom].ele == ele2:
                    continue
                # dislist.append(self.cal_distance(iatom,i))
                dis = self.cal_dis(iatom, i)
                bondtest = bond(self.atom[iatom].ele, self.atom[i].ele, dis)
                bondorder = bondtest.judge_bondorder()
                if bondorder != 0:
                    bondlist.append(dis)
        return bondlist

    def allbondbyele(self, ele1, ele2, flag=1):
        bondall = []
        for iatom in range(self.natom):
            if self.atom[iatom].ele != ele1:
                continue
            bondlist = self.bondsearch(iatom, ele2)
            bondall.extend(bondlist)
        return bondall

    def determinespecies(self, iatom, elelist):
        # dislist = []
        for i in range(self.natom):
            if (i != iatom):
                # dislist.append(self.cal_distance(iatom,i))
                dis = self.cal_dis(iatom, i)
                bondtest = bond(self.atom[iatom].ele, self.atom[i].ele, dis)
                bondorder = bondtest.judge_bondorder()
                if bondorder != 0:  # and self.atom[i].ele == 8:
                    self.atom[iatom].bondlist.append(bondorder)
                    # for id,iza in enumerate(elelist):
                    #    if self.atom[i].ele ==iza:
                    #        if id >0:
                    #            self.atom[iatom].bondtype.append(bondorder+10**id)
                    #        else:
                    #            self.atom[iatom].bondtype.append(bondorder)
                    #        break

                    if self.atom[i].ele == elelist[3]:
                        self.atom[iatom].bondtype.append(bondorder + 1000)
                    if self.atom[i].ele == elelist[2]:
                        self.atom[iatom].bondtype.append(bondorder + 100)
                    if self.atom[i].ele == elelist[1]:
                        self.atom[iatom].bondtype.append(bondorder + 10)
                    if self.atom[i].ele == elelist[0]:
                        self.atom[iatom].bondtype.append(bondorder)
        self.atom[iatom].species = len(self.atom[iatom].bondlist)

        bondall = sum(self.atom[iatom].bondlist)
        self.atom[iatom].bondtype.sort(reverse=True)
        if bondall != 4:
            self.atom[iatom].Ctype = 'radical'
        else:
            if self.atom[iatom].bondtype[-1] == 103:
                self.atom[iatom].Ctype = 'CO'
            elif self.atom[iatom].bondtype[-1] == 102:
                if self.atom[iatom].bondtype[-2] == 101:
                    self.atom[iatom].Ctype = 'acid'
                else:
                    self.atom[iatom].Ctype = 'ketone'
            elif self.atom[iatom].bondtype[-1] == 101:
                self.atom[iatom].Ctype = 'alcohol'
            elif self.atom[iatom].bondtype[-1] == 13:
                self.atom[iatom].Ctype = 'alkyne'
            elif self.atom[iatom].bondtype[-1] == 12:
                self.atom[iatom].Ctype = 'alkene'
            elif self.atom[iatom].bondtype[-1] == 11:
                self.atom[iatom].Ctype = 'alkane'
            elif self.atom[iatom].bondtype[-1] == 1:
                self.atom[iatom].Ctype = 'CH4'

                # self.atom[iatom].species =self.atom[iatom].species+ bondorder
        # dislist.sort()
        # length = min(len(dislist),4)
        # for i in range(length):
        #    if dislist[i]> 1.7:
        #        break
        # self.atom[iatom].species = i + 1
        self.atom[iatom].bondall = bondall
        return

    def calCtypes(self, ):
        self.calOneDimCoord()
        cell = reduce(lambda a, b: a + b, self.lat)
        # print cell
        self.c_natm = pointer(ctypes.c_int(self.natom))
        self.c_xa = pointer(
            (ctypes.c_double * len(self.xa_onedim))(*self.xa_onedim))
        self.c_rv = pointer((ctypes.c_double * len(cell))(*cell))
        iza = [atom.ele for atom in self.atom]
        self.c_iza = pointer((ctypes.c_int * self.natom)(*iza))

    def clear_strC(self):
        self.c_natm = False
        self.c_xa = False
        self.c_rv = False
        self.c_iza = False

    def cal_Q(self):
        self.calCtypes()
        q = ctypes.cdll.LoadLibrary(
            '/home9/shiyf/bin/kpl-path-shiyfgit/Lib_Qcal/calQ.so')

        natom = self.c_natm
        el = pointer(ctypes.c_double(self.energy))
        atom = pointer(ctypes.c_int(0))
        sym = pointer(ctypes.c_int(int(0)))
        za = self.c_iza
        coord = self.c_xa
        rv = self.c_rv
        qglobal = pointer((ctypes.c_double * 4)(*[0.0, 0.0, 0.0, 0.0]))
        # print(natom, za, rv, coord, atom, qglobal, el, sym)
        q.get_order_parameter_(natom, za, rv, coord, atom, qglobal, el, sym)
        qval = list(qglobal.contents)[1:4]
        return [qval[0], qval[1], qval[2]]

    def JudgeBond(self):
        return

    def cal_fragcharge(self):
        "charge here is fakecharge, actually group id"
        groupdict = {}
        for i, atom in enumerate(self.atom):
            # print atom.charge
            try:
                groupdict[atom.charge].append(atom)
            except Exception:
                groupdict[atom.charge] = []
                groupdict[atom.charge].append(atom)
        self.frag = groupdict

    def cal_groupfrag(self):
        groupdict = {}
        for i, atom in enumerate(self.atom):
            # print atom.charge
            try:
                groupdict[self.group[i]].append(atom)
            except Exception:
                groupdict[self.group[i]] = []
                groupdict[self.group[i]].append(atom)
        self.frag = groupdict

    def determinecharge(self):
        self.cal_groupfrag()
        for item in self.frag.values():
            _tmpcharge = 0
            _elelist = []
            for atom in item:
                _elelist.append(atom.ele)
                _tmpcharge = _tmpcharge + atom.ele
            # print _elelist
            if _tmpcharge % 2 != 0:
                for atom in item:
                    if atom.ele == 7:
                        if _elelist == [7, 8]:
                            return 15
                        if _elelist == [7, 8, 8]:
                            return 23
                        return -1
                return 1
        return 2

#    def CheckMin(self,):
#        self.calCtypes()
#        program='/home10/kpl/pymodule/Lib/Lib_fillbond/checkminbond.so'
#        Lminstr = pointer(ctypes.c_bool(0))
#        checkmin = ctypes.cdll.LoadLibrary(program)
#        checkmin.judgebond_(self.c_natm,self.c_iza,self.c_xa,self.c_rv,Lminstr)
#        self.Lminstr = bool(Lminstr.contents)
#        return bool(Lminstr.contents)

    def CheckMin(self, flag=1):
        sqnatm = self.natom**2
        self.calCtypes()
        # print("use Lib_fillbond_tmp")
        # program = '/home10/kpl/pymodule/Lib/Lib_fillbond/checkminbond.so'
        program = "/home9/shiyf/bin/kpl-path-shiyfgit/Lib_fillbond_tmp/checkminbond.so"
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
            checkmin.judgebondsurface_(self.c_natm, self.c_iza, self.c_xa,
                                       self.c_rv, Lminstr, bmatrix, bondneed,
                                       surface)

        bmx = list(bmatrix.contents)
        # self.bmx2D = np.array(bmx).reshape(self.natom, self.natom)
        # self.bmx1D = bmx

        self.bondneed = list(bondneed.contents)
        self.Lminstr = bool(Lminstr.contents)

        return self.Lminstr, bmx, list(bondneed.contents), list(
            surface.contents)

    def Bondmatrix(self, ):

        program = '/home10/kpl/pymodule/Lib/Lib_bondmatrix/bondmatrix.so'
        "/home9/shiyf/bin/kpl-path-shiyfgit/Lib_fillbond_tmp/bondmatrix.so"
        self.calCtypes()
        sqnatm = self.natom**2
        self.cdnt2fcnt()
        onedim_fa = reduce(lambda a, b: a + b, self.fdnt)
        self.c_fxa = pointer((ctypes.c_double * len(onedim_fa))(*onedim_fa))

        bmatrix = pointer((ctypes.c_int * sqnatm)(*[0 for i in range(sqnatm)]))
        bcal = ctypes.cdll.LoadLibrary(program)
        bcal.hbondmatrix_(self.c_natm, self.c_fxa, self.c_iza, self.c_rv,
                          bmatrix)
        return list(bmatrix.contents)

    def Segmolecular(self, recal=True):
        """ warning : bond matrix must first be calculated."""
        program = '/home10/kpl/pymodule/Lib/Lib_bondmatrix/bondmatrix.so'
        if recal:
            self.calCtypes()
            self.cdnt2fcnt()
            onedim_fa = reduce(lambda a, b: a + b, self.fdnt)
            self.c_fxa = pointer(
                (ctypes.c_double * len(onedim_fa))(*onedim_fa))
        self.c_bmx = pointer((ctypes.c_int * len(self.bmx1D))(*self.bmx1D))
        group = pointer(
            (ctypes.c_int * self.natom)(*[0 for i in range(self.natom)]))
        scal = ctypes.cdll.LoadLibrary(program)
        scal.hsegmentmol_(self.c_natm, self.c_fxa, self.c_iza, self.c_rv,
                          group, self.c_bmx)
        # print(list(group.contents))
        return list(group.contents)

    def GetPatternAtom(self, pattern):
        Atom1 = pattern[0]
        Atom2 = pattern[1]
        mode = int(pattern[-1])

        atomlist1 = []
        design = re.findall(r"\d+\.?\d*", Atom1)
        if len(design) != 0:
            for item in design:
                atomlist1.append(int(item) - 1)
        else:
            for iatom, atom in enumerate(self.atom):
                if atom.elesymbol == Atom1:
                    atomlist1.append(iatom)

        atomlist2 = []
        design = re.findall(r"\d+\.?\d*", Atom2)
        if len(design) != 0:
            for item in design:
                atomlist2.append(int(item) - 1)
        else:
            atomlist2 = []
            for iatom, atom in enumerate(self.atom):
                if atom.elesymbol == Atom2:
                    atomlist2.append(iatom)

        return atomlist1, atomlist2, mode

    def GetAllBond(self):
        self.allbond = []
        for i in range(self.natom):
            for j in range(i + 1, self.natom):
                if self.bmx2D[i][j] > 0:
                    self.allbond.append([i, j, self.bmx2D[i][j]])

        self.nbond = len(self.allbond)

    def GetAtomInfo(self):
        for i in range(self.natom):
            self.atom[i].expbond = 0
            self.atom[i].imph = 0
        # print(self.bmx2D)
        for i in range(self.natom):
            for j in range(i + 1, self.natom):
                if self.bmx2D[i][j] > 0:
                    if self.atom[j].ele != 1:
                        self.atom[i].expbond = self.atom[i].expbond + 1
                    else:
                        self.atom[i].imph = self.atom[i].imph + 1
                    if self.atom[i].ele != 1:
                        self.atom[j].expbond = self.atom[j].expbond + 1
                    else:
                        self.atom[j].imph = self.atom[j].imph + 1

    def CalUnsaturatedNumber(self):
        nH = 0
        needH = 2
        nheavyatom = 0
        for atom in self.atom:
            if atom.ele == 1:
                nH = nH + 1
            else:
                nheavyatom = nheavyatom + 1
                needH = needH + (8 - atom.ele)
        self.UnsNum = (needH - nH) / 2
        self.nheavyatom = nheavyatom

    def RemoveMetalBond(self):
        self.surfaceatom = [0 for i in range(self.natom)]
        self.bondneed = [0 for i in range(self.natom)]
        self.bmxsave = np.zeros((self.natom, self.natom), dtype=np.int)

        for i in range(self.natom):
            for j in range(self.natom):
                if self.bmx2D[i][j] > 0:
                    self.bmxsave[i][j] = self.bmx2D[i][j]
                    if self.atom[i].ele > 18:  # or self.atom[j].ele > 18:
                        self.bmx2D[i][j] = 0
                        self.surfaceatom[j] = 1
                    if self.atom[j].ele > 18:
                        self.bmx2D[i][j] = 0
                        self.surfaceatom[i] = 1

        self.bmx1D = []
        for line in self.bmx2D:
            self.bmx1D.extend(line)

    def ChangeMetalAtom(self, elesymbol1, elesymbol2, expand=1):
        ele1 = PT.Eledict[elesymbol1]
        ele2 = PT.Eledict[elesymbol2]
        # newstr = Str()

        self.cdnt2fcnt()
        self.abc[0] = self.abc[0] * expand
        self.abc[1] = self.abc[1] * expand
        self.abc2lat()

        for iatom, atom in enumerate(self.atom):
            if atom.ele == ele1:
                atom.ele = ele2
                atom.elesymbol = elesymbol2
                atom.xyz = np.matmul(self.fdnt[iatom], self.lat)

        return

    def GetAllNeighbours(self, rcutoff):
        self.cal_dis_matrix()
        self.allnbr = []
        for i in range(self.natom):
            nbr = []
            for j in range(self.natom):
                if self.Dmatrix[i][j] < rcutoff and i != j:
                    nbr.append([i, self.Dmatrix[i][j], j])
            self.allnbr.append(nbr)

        return self.allnbr

    def PeriodicInfo(self):
        self.cdnt2fcnt()
        fdnt = self.fdnt[:]
        for i, x in enumerate(fdnt):
            fdnt[i] = map(lambda y: (y + 1000.0) % 1.0, x)

        cdnt = np.dot(fdnt, self.lat)

        shortdis = {}

        for i in range(0, self.natom):
            cart0 = cdnt[i]
            for j in range(i + 1, self.natom):
                b = np.subtract(cdnt[j], cart0)
                d = m.sqrt(np.dot(b, b))
                shortdis[i, j] = [d, [0, 0, 0]]
                shortdis[j, i] = [d, [0, 0, 0]]

        cut1 = [-1, 0, 1]
        for i0 in cut1:
            for i1 in cut1:
                for i2 in cut1:
                    if i0 == 0 and i1 == 0 and i2 == 0:
                        continue
                    # print i0,i1,i2
                    frac2 = []
                    for i in range(self.natom):
                        frac2.append(np.add(fdnt[i], [i0, i1, i2]))
                    # print frac2
                    cart2 = np.dot(frac2, self.lat)
                    for i in range(0, self.natom):
                        cart0 = cdnt[i]
                        for j in range(i + 1, self.natom):
                            b = np.subtract(cart2[j], cart0)
                            d = m.sqrt(np.dot(b, b))
                            if d < shortdis[i, j][0]:
                                shortdis[i, j] = [d, [i0, i1, i2]]
                                shortdis[j, i] = [d, [-i0, -i1, -i2]]
        self.supercellinfo = shortdis

#    def GetNeighbour(self,iatom):
#        for i in range(self.natom):
#            if self.bmx2D[iatom][i]>0:
#                self.atom[iatom].neighbour.append()
#
#    def GetAtomnodes(self):
#
#        for iatom,atom in enumerate(self):
#            self[i]

#
# "=============== the following part is designed for transfer ============="
# "must excute transfer before calling the different part function"
#

    def TransferToXYZcoordStr(self):
        self.Nat = self.natom
        self.Latt = self.abc
        self.Cell = self.Latt2Cell()
        self.Energy = self.energy

        self.EleNam = []
        self.iza = []
        self.Coord = []
        for atom in self.atom:
            self.EleNam.append(atom.elesymbol)
            self.iza.append(atom.ele)
            self.Coord.append(atom.xyz)

        for iele, elesymbol in enumerate(self.elenameList):
            self.sp[elesymbol] = self.natompe[iele]
            self.sporder[elesymbol] = iele + 1

        if self.Lfor:
            self.For = []
            for atom in self.atom:
                self.For.append(atom.force)
            self.GetMaxF()
            self.maxF = self.maxf

    def TransferToKplstr(self):
        self.atom = []
        for i in range(self.Nat):
            self.atom.append(S_atom(self.Coord[i], self.iza[i]))
            if self.Lfor:
                self.atom[i].force = self.For[i]
        self.abc = self.Latt
        self.energy = self.Energy
        self.sortatombyele()
        self.calAtomnum()
        self.abc2lat()

#
# "=============== the following part is copyed from zpliu XYZCoord.py ============="
# "must excute TransferToXYZcoordStr before calling the following function"
#

    def printCell(self):
        """ formated cell print---used to generate INPUT_DEBUG"""
        s = ""
        for x in self.Cell:
            s += "%10.4f %10.4f %10.4f\n  " % (x[0], x[1], x[2])
        return s

    def printCoord(self):
        """ used to generate INPUT_DEBUG"""
        s = ""
        for i, xa in enumerate(self.Coord):
            s += " %15.9f  %15.9f  %15.9f   %3d\n" % (
                xa[0], xa[1], xa[2], self.sporder[self.EleNam[i]])
        return s

    def myfilter(self, BadStrP):
        if (self.Energy > BadStrP.HighE or self.Energy < BadStrP.LowE
                or min(self.Latt[:3]) < BadStrP.MinLat
                or max(self.Latt[:3]) > BadStrP.MaxLat
                or max(self.Latt[3:7]) > BadStrP.MaxAngle
                or min(self.Latt[3:7]) < BadStrP.MinAngle
                #          max(self.Latt[:3]) > 5.0*min(self.Latt[:3]) or \
                or self.maxFF > BadStrP.MaxFor):
            # print '111',self.maxFF,BadStrP.MaxFor
            return False
        else:
            # print self.Energy,True
            # return True
            if len(self.For):
                # print self.maxF,BadStrP.MaxFor
                if self.maxF > BadStrP.MaxFor:
                    return False
                else:
                    return True
            else:
                return True

    def myfilter_byd(self, Badd1, Badd2, type):

        L = True
        x1 = [x for x in Badd1.keys()]
        y1 = [y for y in Badd2.keys()]
        #        print 'XXX',x1, y1
        for key, value in self.d.items():
            #           print key, value, Badd1[key], Badd2[key]
            if key in x1:
                if type == 1 and value < Badd1[key]:
                    L = False
                    continue
            if key in y1:
                if type == 1 and value > Badd2[key]:
                    L = False
                    continue
            if key in x1 and key in y1:
                if type == 2 and value > Badd1[key] and value < Badd2[key]:
                    L = False
        return L

    def Latt2Cell(self):
        """ transform of (a,b,c,alpha,beta,gamma) to lattice vector"""
        # lower-triangle
        a, b, c, alpha, beta, gamma = self.Latt[:6]
        pi = 3.14159265
        alpha, beta, gamma = alpha * pi / 180.0, beta * pi / 180.0, gamma * pi / 180.0
        bc2 = b**2 + c**2 - 2 * b * c * m.cos(alpha)
        h1 = a
        h2 = b * m.cos(gamma)
        h3 = b * m.sin(gamma)
        h4 = c * m.cos(beta)
        h5 = ((h2 - h4)**2 + h3**2 + c**2 - h4**2 - bc2) / (2 * h3)
        # x = c**2 - h4**2 - h5**2
        h6 = m.sqrt(c**2 - h4**2 - h5**2)
        cell = [[h1, 0., 0.], [h2, h3, 0.], [h4, h5, h6]]
        return cell

    def Cell2Latt(self):
        """ transform of lattice vector to (a,b,c,alpha,beta,gamma)"""
        lat = self.Cell[0:3]
        a = np.linalg.norm(lat[0])
        b = np.linalg.norm(lat[1])
        c = np.linalg.norm(lat[2])
        pi = 3.14159265
        alpha = m.acos(np.dot(lat[1], lat[2]) / (b * c)) * 180.0 / pi
        beta = m.acos(np.dot(lat[0], lat[2]) / (a * c)) * 180.0 / pi
        gamma = m.acos(np.dot(lat[0], lat[1]) / (a * b)) * 180.0 / pi
        return [a, b, c, alpha, beta, gamma]

    def Volume(self):
        """ calculate volume of structure"""
        a = self.Cell[0]
        b = self.Cell[1]
        c = self.Cell[2]
        return np.dot(np.cross(a, b), c)

    def ReciCell(self):
        """ reciprocal lattice """
        return np.transpose(np.linalg.inv(self.Cell))

    def FracCoord(self):
        """ fractional coordinate """
        cellr = np.linalg.inv(self.Cell)
        return np.dot(self.Coord, cellr)

    def frac_modulo(self):
        frac = self.FracCoord()
        for i, x in enumerate(frac):
            frac[i] = map(lambda y: (y + 1000.0) % 1.0, frac[i])
        return frac

    def centralize(self, n=1):
        if n == 1:
            frac = self.frac_modulo()
        else:
            #           key= self.centerf.keys() ; key.sort()       #  sort by keys in dict
            #           frac = np.array([self.centerf[i] for i in key])
            frac = [self.centerf[i] for i in sorted(self.centerf.keys())]
#       print frac
        cart = np.dot(frac, self.Cell)

        # self.Cell=[[10,0,0],[0,10,0],[0,0,10]]

        self.pos = np.dot([0.5, 0.5, 0.5], self.Cell)
        com = []
        for i in range(3):
            com.append(sum([x[i] for x in cart]) / float(self.Nat))
        self.cart = np.array(
            [np.add(np.subtract(x, com), self.pos) for x in cart])
        self.frac = np.dot(self.cart, np.linalg.inv(self.Cell))
#       return self.frac

    def frac_center(self, at):
        # recursive to position all fractional coordinate
        if len(self._cycle) == 0:
            self.centerf[at.keys()[0]] = at.values()[0]
            self._currentlist = []
            self._currentlist.append(at.keys()[0])
        self._cycle.append(at.keys()[0])

        if len(self.centerf) == self.Nat:
            return True
        else:
            neig = self.neighbor(at)
            # get frac of atom at  neig ={1:[XX,XX,XX],2:[XX,XX,XX]}
            for i in neig.keys():
                if i not in self.centerf.keys():
                    self.centerf[i] = neig[i]
                    self._currentlist.append(i)
                    if len(self.centerf) == self.Nat:
                        return True
            for i in neig.keys():
                if i not in self._cycle:
                    self.frac_center({i: self.centerf[i]})
                if len(self.centerf) == self.Nat:
                    return True
            if at.keys()[0] == 0:
                return False

    def ChemicalFormula(self):
        # recursive to position all fractional coordinate
        Nmol = 0
        fragment = {}
        while len(self.centerf) != self.Nat:
            for i in range(self.Nat):
                if i not in self.centerf.keys():
                    self._cycle = []
                    self._currentlist = []
                    L = self.frac_center({i: self.FracCoord()[i]})
                    Nmol += 1
                    elelist = {}
                    fragment[Nmol] = ""
                    # print self._currentlist, len(self._currentlist)
                    for j in self._currentlist:
                        if self.EleNam[j] in elelist.keys():
                            elelist[self.EleNam[j]] += 1
                        else:
                            elelist[self.EleNam[j]] = 1
                    for j in elelist.keys():
                        fragment[Nmol] = fragment[Nmol] + j + str(elelist[j])
                    if L:
                        self.centralize(0)
                        for j in range(self.Nat):
                            self.Coord[j][0:3] = self.cart[j][0:3]
                        return fragment, self.Coord, self.Cell

    def JudgeShape(self, cut=2.6):
        """ find the shortest bond between neighboring cell to judge solid shape"""
        if len(self.centerf) == 0:
            L = self.frac_center({0: self.FracCoord()[0]})
        if not L:
            return {"fragments": -1}

        # vect = {}
        # vect_1 = {}
        com = []
        cellconnect = {}
        self.centralize(0)

        V = []
        # vect = {}
        Vcom = {}
        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    if (i == 0 and j == 0 and k == 0):
                        continue
                    V.append([i, j, k])
                    com = np.dot(np.add([0.5, 0.5, 0.5], [i, j, k]), self.Cell)
                    Vcom[str(i) + str(j) + str(k)] = np.subtract(com, self.pos)

        line = False
        layer = False
        # solid = False
        # count = 0
        normalvector = []
        cellconnect = {}
        for i0, i1, i2 in V:
            if not str(i0) + str(i1) + str(i2) in cellconnect.keys():
                cellconnect[str(i0) + str(i1) +
                            str(i2)] = self.cell_bondconnect(
                                cut, [i0, i1, i2])
            if cellconnect[str(i0) + str(i1) + str(i2)]:
                line = True
                layer_hole_detect = []
                for j0, j1, j2 in V:
                    if i0 == j0 and i1 == j1 and i2 == j2:
                        continue
                    angle = self.vect_angle(Vcom[str(i0) + str(i1) + str(i2)],
                                            Vcom[str(j0) + str(j1) + str(j2)])
                    if abs(angle - 90) < 70:
                        if not str(j0) + str(j1) + str(j2) in cellconnect.keys(
                        ):
                            cellconnect[str(j0) + str(j1) +
                                        str(j2)] = self.cell_bondconnect(
                                            cut, [j0, j1, j2])
                        if cellconnect[str(j0) + str(j1) + str(j2)]:
                            layer_hole_detect.append(angle)
                            layer = True
                            VectorC = np.cross(
                                Vcom[str(i0) + str(i1) + str(i2)],
                                Vcom[str(j0) + str(j1) + str(j2)])
                            for i in range(len(normalvector) - 1):
                                for j in range(i + 1, len(normalvector)):
                                    angle = self.vect_angle(
                                        normalvector[i], normalvector[j])
                                    if angle > 40:
                                        return {
                                            "solid": 7777
                                        }  # early break if two intercected layer identified
                            angle = 0.0
                            for i in range(len(normalvector)):
                                angle = max(
                                    angle,
                                    self.vect_angle(normalvector[i], VectorC))
                            if len(normalvector) > 1 and angle < 10:
                                break
                            normalvector.append(VectorC)
                            hole_detect = []
                            for k0, k1, k2 in V:
                                if k0 == i0 and k1 == i1 and k2 == i2:
                                    continue
                                if k0 == j0 and k1 == j1 and k2 == j2:
                                    continue
                                angle = self.vect_angle(
                                    Vcom[str(k0) + str(k1) + str(k2)], VectorC)
                                if angle < 70:
                                    if not str(k0) + str(k1) + str(
                                            k2) in cellconnect.keys():
                                        cellconnect[
                                            str(k0) + str(k1) +
                                            str(k2)] = self.cell_bondconnect(
                                                cut, [k0, k1, k2])
                                    if cellconnect[str(k0) + str(k1) +
                                                   str(k2)]:
                                        hole_detect.append(angle)
                            # ------THIRD LOOP END----
                            # print hole_detect
                            if len(hole_detect) == 0:
                                return {"layer": 7777}
                            elif (max(hole_detect) - min(hole_detect)) < 15:
                                return {"solid-largehole": 9999}
                            elif (max(hole_detect) - min(hole_detect)) < 35:
                                return {"solid-smallhole": 8888}
                            else:
                                return {"solid-dense": 6666}

                #        else: continue  ------SECOND LOOP END----
                # the following info on layer should not be met due to the early return above
                if len(layer_hole_detect) == 0:
                    return {"line": 7777}
                elif (max(layer_hole_detect) - min(layer_hole_detect)) < 15:
                    return {"layer-largehole": 9999}
                elif (max(layer_hole_detect) - min(layer_hole_detect)) < 35:
                    return {"layer-smallhole": 8888}
                else:
                    return {"layer-dense": 6666}
        #   else: continue   ---- FIRST LOOP END------

        if layer:
            return {"?layer?": 6666}  # should not be useful
        elif line:
            return {"?line?": 6666}  # should not be useful
        else:
            return {"cluster": 6666}

    def vect_angle(self, V0, V1):
        a = np.dot(V0, V1) / m.sqrt(np.dot(V0, V0)) / m.sqrt(np.dot(V1, V1))
        angle = 180.0
        if abs(abs(a) - 1.0) > 1e-6:
            angle = m.acos(a) * 180.0 / 3.14159265
        if angle > 90.0:
            angle = 180.0 - angle
        return angle

    def cell_bondconnect(self, cut, cell1, cell0=[0, 0, 0]):
        frac1 = []
        # VectorA = []
        N = self.Nat
        for i in range(N):
            frac1.append(np.add(self.frac[i], cell1))
        cart1 = np.dot(frac1, self.Cell)

        if cell0[0] != 0 or cell0[1] != 0 or cell0[2] != 0:
            frac0 = []
            for i in range(N):
                frac0.append(np.add(self.frac[i], cell0))
            cart0 = np.dot(frac0, self.Cell)
            pos = []
            pos.append(sum([y[0] for y in cart0]) / float(self.Nat))
            pos.append(sum([y[1] for y in cart0]) / float(self.Nat))
            pos.append(sum([y[2] for y in cart0]) / float(self.Nat))
        else:
            # fract0 = self.frac
            cart0 = self.cart
            pos = self.pos

        d0 = 999
        for i in range(N):
            for j in range(N):
                b = np.subtract(cart1[i], cart0[j])
                d = m.sqrt(np.dot(b, b))
                if d < d0:
                    d0 = d
        return d0 < cut

    def neighbor(self, at, cut=2.6):
        """ find the shortest bond neighbors of at  """
        frac = self.frac_modulo()
        cart = np.dot(frac, self.Cell)
        cart0 = np.dot(at.values()[0], self.Cell)
        N = self.Nat  # len(cart)
        neig = {}
        neig2 = {}
        for j in range(N):
            if j == at.keys()[0]:
                continue
            b = np.subtract(cart0, cart[j])
            d = m.sqrt(np.dot(b, b))
            if d < cut:
                neig[j] = frac[j]
                neig2[j] = d
                # print d,neig
            # print bonddict
        cut1 = [-1, 0, 1]
        for i0 in cut1:
            for i1 in cut1:
                for i2 in cut1:
                    if i0 == 0 and i1 == 0 and i2 == 0:
                        continue
                    # print i0,i1,i2
                    frac2 = []
                    for i in range(N):
                        frac2.append(np.add(frac[i], [i0, i1, i2]))
                    cart2 = np.dot(frac2, self.Cell)
                    for j in range(N):
                        b = np.subtract(cart2[j], cart0)
                        d = m.sqrt(np.dot(b, b))
                        if d < cut:
                            if j in neig.keys():
                                if neig2[j] > d:
                                    neig[j] = frac2[j]
                                    neig2[j] = d
                            else:
                                neig[j] = frac2[j]
                                neig2[j] = d

                            # print d,neig
        # neig = sorted(neig,keys=neigh.valueskeys())
        return neig

    def Shortestbond(self):
        """ find the shortest bond """
        frac = self.frac_modulo()
        cart = np.dot(frac, self.Cell)
        N = self.Nat  # len(cart)
        bonddict = {}
        if N == 1:
            bondname = self.EleNam[0] + "-" + self.EleNam[0]
            bonddict[bondname] = 10.0
            return bonddict
        for i in range(N):
            for j in range(N):
                bondname = self.EleNam[i] + "-" + self.EleNam[j]
                #               if self.sporder[self.EleNam[i]]>self.sporder[self.EleNam[j]]: bondname=self.EleNam[j]+"-"+self.EleNam[i]
                if PT.Eledict[self.EleNam[i]] > PT.Eledict[self.EleNam[j]]:
                    bondname = self.EleNam[j] + "-" + self.EleNam[i]
                bonddict[bondname] = 10.0
                # return bonddict
        for i in range(N):
            for j in range(i + 1, N):
                b = np.subtract(cart[i], cart[j])
                d = m.sqrt(np.dot(b,
                                  b))  # (cart[i]-cart[j]),(cart[i]-cart[j])))
                bondname = self.EleNam[i] + "-" + self.EleNam[j]
                # print bondname
                #               if self.sporder[self.EleNam[i]]>self.sporder[self.EleNam[j]]: bondname=self.EleNam[j]+"-"+self.EleNam[i]
                if PT.Eledict[self.EleNam[i]] > PT.Eledict[self.EleNam[j]]:
                    bondname = self.EleNam[j] + "-" + self.EleNam[i]
                if bondname not in bonddict:
                    bonddict[bondname] = d
                else:
                    bonddict[bondname] = min(d, bonddict[bondname])
                # print bonddict


# for cluster
#       return bonddict
# for cluster
        cut1 = [-1, 0, 1]
        for i0 in cut1:
            for i1 in cut1:
                for i2 in cut1:
                    if i0 == 0 and i1 == 0 and i2 == 0:
                        continue
                    # print i0,i1,i2
                    frac2 = []
                    for i in range(N):
                        frac2.append(np.add(frac[i], [i0, i1, i2]))
                    cart2 = np.dot(frac2, self.Cell)
                    for i in range(N):
                        for j in range(N):
                            b = np.subtract(cart2[i], cart[j])
                            d = m.sqrt(np.dot(
                                b,
                                b))  # (cart2[i]-cart[j]),(cart2[i]-cart[j])))
                            bondname = self.EleNam[i] + "-" + self.EleNam[j]
                            # if self.sporder[self.EleNam[i]]>self.sporder[self.EleNam[j]]: bondname=self.EleNam[j]+"-"+self.EleNam[i]
                            if PT.Eledict[self.EleNam[i]] > PT.Eledict[
                                    self.EleNam[j]]:
                                bondname = self.EleNam[j] + "-" + self.EleNam[i]
                            if bondname not in bonddict:
                                bonddict[bondname] = d
                            bonddict[bondname] = min(d, bonddict[bondname])
                            # print bonddict['Ge-Ge']
        return bonddict

    def _Gen_arc(self, coord, fname='outtmp.arc'):
        with open(fname, 'w') as fout:
            fout.write("!BIOSYM archive 2\nPBC=ON\n")
            energy = self.Energy
            i = 0
            if not self.Lfor:
                fout.write("\t\t\t\tEnergy\t%8d        -0.0000  %17.6f\n" %
                           (i + 1, energy))
            else:
                fout.write("\t\t\t\tEnergy\t%8d        %10.6f  %17.6f\n" %
                           (i + 1, self[i].maxF, energy))

            fout.write("!DATE\n")
            lat = self.Latt[0:6]
            fout.write("PBC  %12.6f%12.6f%12.6f%12.6f%12.6f%12.6f\n" %
                       (lat[0], lat[1], lat[2], lat[3], lat[4], lat[5]))
            for j in range(self.Nat):
                ele = self.EleNam[j]
                xa = coord[j]
                fout.write(
                    "%-2s%18.9f%15.9f%15.9f CORE %4d %-2s %-2s   0.0000 %4d\n"
                    % (ele, xa[0], xa[1], xa[2], j + 1, ele, ele, j + 1))
            fout.write("end\nend\n")

    def hkl_dspacking(self, hkl):

        hklLib = ctypes.cdll.LoadLibrary('./crystal-hkl-subroutines.so')
        onedimrv = reduce(lambda a, b: a + b, self.Cell)
        rv = pointer((ctypes.c_double * len(onedimrv))(*onedimrv))
        hkl0 = pointer((ctypes.c_int * len(hkl))(*hkl))
        d = 0.00
        d_c = pointer((ctypes.c_double(d)))

        hklLib.plane_dspacing2_(rv, hkl0, d_c)

        d = float(str(d_c.contents).split('(')[1].split(')')[0])

        return d

    def SteinhartQ_cal(self):

        q = ctypes.cdll.LoadLibrary('./calQ.so')
        xa = reduce(lambda a, b: a + b, self.Coord)
        cell = reduce(lambda a, b: a + b, self.Cell)

        natom = pointer(ctypes.c_int(self.Nat))
        el = pointer(ctypes.c_double(self.Energy))
        atom = pointer(ctypes.c_int(0))
        sym = pointer(ctypes.c_int(int(0)))
        za = pointer((ctypes.c_int * len(self.iza))(*self.iza))
        coord = pointer((ctypes.c_double * len(xa))(*xa))
        rv = pointer((ctypes.c_double * len(cell))(*cell))
        qglobal = pointer((ctypes.c_double * 4)(*[0.0, 0.0, 0.0, 0.0]))
        q.get_order_parameter_(natom, za, rv, coord, atom, qglobal, el, sym)
        qval = (qglobal.contents)[1:4]
        return (qval[0], qval[1], qval[2])


class FooError(Exception):
    pass


class BadStr(object):
    def __init__(self):
        self.HighE = -0
        self.LowE = -9999999
        self.MinAngle = 50
        self.MaxAngle = 130
        self.MinLat = 1.6
        self.MaxLat = 50.0
        self.MaxFor = 10.0


def ParaWrap_JudgeShape(x):
    return x.JudgeShape()


def ParaWrap_Shortestbond(x):
    return x.Shortestbond()


def ParaWrap_ChemicalFormula(x):
    return x.ChemicalFormula()


def ParaWrap_SteinhartQ_cal(x):
    return x.SteinhartQ_cal()


if __name__ == '__main__':
    a = 1
    print(a)
