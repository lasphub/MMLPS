from ECFP import allstr
from structure_new import Str
import os
import numpy as np
# from functools import reduce


class AllStr(allstr):
    def splituncm(self):
        workdirs = []
        for i in range(0, len(self), 2):
            # if self[i].upper ==1 or self[i+1].upper ==1:
            #    print  ('upper surface wrong path-%d'%(i/2))
            #    continue
            workdir = 'path-%d' % (i / 2)
            workdirs.append(workdir)
            os.mkdir(workdir)
            reactpair = allstr()
            reactpair.append(self[i])
            reactpair.printall('%s/input.arc' % (workdir))
            reactpair.append(self[i + 1])
            reactpair.printall('%s/uncm.arc' % (workdir))
            os.system('cp ./input %s' % (workdir))
            os.system('ln -s %s/*.pot %s' % (os.getcwd(), workdir))
        return workdirs

    def singlesplit(self):
        workdirs = []
        for i in range(0, len(self)):
            workdir = 'para%d' % (i)
            workdirs.append(workdir)
            os.mkdir(workdir)
            reactpair = allstr()
            reactpair.append(self[i])
            reactpair.printall('%s/input.arc' % (workdir))
            os.system('cp ./input %s' % (workdir))
            os.system('ln -s %s/*.pot %s' % (os.getcwd(), workdir))
        return workdirs

    def readfile(self, inputfile, forcefile=False):
        f = open(inputfile, 'r')
        currentStr = -1
        for line in f:
            if ('Energy' in line or 'React' in line
                    or ('TS' in line and 'sect' not in line) or 'Str' in line):
                self.append(Str())
                currentStr += 1
                try:
                    self[currentStr].energy = float(line.split()[-1])
                except BaseException:
                    self[currentStr].energy = float(line.split()[-2])
            elif 'CORE' in line:
                self[currentStr].addatom(line, 1)
            elif ('PBC' in line) and ('ON' not in line):
                self[currentStr].abc = [float(x) for x in line.split()[1:]]
        f.close()
        for str in self:
            str.calAtomnum()
            str.abc2lat()
            # str.screenuppersurf()

        # self.GetAllsminame()
        # self.calAllFakeBondMaxtrix()

        if forcefile:
            f = open(forcefile, 'r')
            currentStr = -1
            for line in f:
                if "For" in line:
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

    def RPpairmode(self):
        return

    def GenNameDict(self, filename):
        self.GetAllName(filename)
        self.namedict = {}
        wrongstr = []
        for i in range(len(self)):
            # if self[i].upper ==1: continue
            if not self[i].Lminstr:
                wrongstr.append(i)
                continue
            name = self[i].sminame
            if name not in self.namedict.keys():
                self.namedict[name] = [i]
            else:
                self.namedict[name].append(i)
        self.nstrpername = [(x, len(y)) for x, y in self.namedict.items()]
        self.nstrpername.sort(key=lambda X: X[1], reverse=True)
        self.printlist(wrongstr, 'wrongmin.arc')
        os.system('cat wrongmin.arc >> ../wrongmin.arc')

    def ScreenPattern(self, pattern, inistr=0):
        atomlist1, atomlist2, mode = self[0].GetPatternAtom(pattern)
        bmxi = np.array(self[inistr].Bondmatrix()).reshape(
            self[inistr].natom, self[inistr].natom)
        pproduct = []
        for name in self.namedict.keys():
            istr = self.namedict[name][0]
            bmxj = np.array(self[istr].Bondmatrix()).reshape(
                self[istr].natom, self[istr].natom)
            lpattern = self.JudgePattern(bmxi, bmxj, atomlist1, atomlist2,
                                         mode)
            if lpattern:
                pproduct.append(name)
                print(name)
        return pproduct

    def JudgePattern(self, bmxi, bmxj, atomlist1, atomlist2, mode):
        lpattern = 0
        change = bmxj - bmxi
        # atomlist1,atomlist2,mode=self[0].GetPatternAtom(pattern)
        for i in atomlist1:
            for j in atomlist2:
                if change[i][j] == mode:
                    lpattern = 1
                    break
        return lpattern

    def SelectStr(self, pattern, trace):
        nopattern = 0
        self.GenNameDict('uncm.arc')
        pproduct = self.ScreenPattern(pattern)
        if len(pproduct) != 0 and nopattern == 0:
            pplist = [(name, self.namedict[name]) for name in pproduct]
        else:
            print('No match pattern products')
            pplist = [(name, self.namedict[name])
                      for name in self.namedict.keys()]
        pplist.sort(key=lambda X: len(X[1]), reverse=True)
        for name, namelist in pplist:
            if name not in trace:
                outlist = namelist
                trace.append(name)
                break
        # print(outlist)
        out = AllStr()
        for i in outlist:
            out.append(self[i])
        out.sort(key=lambda x: x.energy)
        out.printlist([0], 'next.arc')
        return trace


#    def GetReactCenter(self, depth =1):
#        # must have caled bmx
#
#        str1 = self[0]
#        str2 = self[1]
#
#        str1.GetAtomInfo()
#        str2.GetAtomInfo()
#
#        bmxc = str1.bmx2D - str2.bmx2D
#        atomlist = []
#        for i in range(str1.natom):
#            for j in xrange(i,str1.natom):
#                if bmxc[i][j] != 0:
#                    if i not in atomlist:  atomlist.append(i)
#                    if j not in atomlist:  atomlist.append(j)
#
#        atomlist0= atomlist[:]
#        if depth ==1:
#            for i in atomlist0:
#                for j in range(str1.natom):
#                    if str1.bmx2D[i][j] >0  and j not in atomlist :  atomlist.append(j)
#                    if str2.bmx2D[i][j] >0  and j not in atomlist :  atomlist.append(j)
#
#        #atomlist.sort()
#
#        _tmpstr1 = Str()
#        _tmpstr2 = Str()
#        rcatomindice= {}
#
#        for i in atomlist:
#            if str1.atom[i].ele == 1: continue
#            _tmpstr1.atom.append(str1.atom[i])
#            _tmpstr2.atom.append(str2.atom[i])
#            rcatomindice[len(_tmpstr1.atom)-1]=i
#
#        _tmpstr1.abc = str1.abc[:]
#        _tmpstr2.abc = str2.abc[:]
#
#        _tmpstr1.calAtomnum()
#        _tmpstr1.abc2lat()
#        _tmpstr2.calAtomnum()
#        _tmpstr2.abc2lat()
#
#        _tmpstr1.bmx2D=np.array([[0 for i in range(_tmpstr1.natom)] for i in range(_tmpstr1.natom)])
#        _tmpstr2.bmx2D=np.array([[0 for i in range(_tmpstr1.natom)] for i in range(_tmpstr1.natom)])
#
#        for i in range(_tmpstr1.natom):
#            for j in range(_tmpstr1.natom):
#                _tmpstr1.bmx2D[i][j]= str1.bmx2D[rcatomindice[i]][rcatomindice[j]]
#                _tmpstr2.bmx2D[i][j]= str2.bmx2D[rcatomindice[i]][rcatomindice[j]]
#
#
#        for i,atom in enumerate(_tmpstr1.atom):
#            atom.expbond = str1.atom[rcatomindice[i]].expbond
#            atom.imph =  str1.atom[rcatomindice[i]].imph
#            atom.group = str1.group[i]
#
#        for i,atom in enumerate(_tmpstr2.atom):
#            atom.expbond = str2.atom[rcatomindice[i]].expbond
#            atom.imph =  str2.atom[rcatomindice[i]].imph
#            atom.group = str2.group[i]
#
#        reactcenter = allstr()
#        reactcenter.append(_tmpstr1)
#        reactcenter.append(_tmpstr2)
#
#        reactcenter.OutMolFile(printlist= [0], outfile = 'rc1.mol')
#        reactcenter.OutMolFile(printlist= [1], outfile = 'rc2.mol')
#
#        reactcenter.printall('reactcenter.arc')
#
#        return reactcenter
#
