#! /usr/bin/env python

import sys
import os
# import numpy
# from atom_k import S_atom
# import time
from structure_new import Str
import multiprocessing
from multiprocessing import Pool
# import glob
# import shutil
import traceback
import subprocess
import ctypes
from ctypes import pointer
from functools import reduce


class Hostfile(object):
    def __init__(self, workdir, cpuperjob, setmasternode=0):
        self.workdir = workdir
        self.cpuperjob = cpuperjob
        self.setmasternode = setmasternode

    def setHostfile(self):
        hostDict = self.getProc()

        divHost, totalproc = self.alloProc(hostDict, self.cpuperjob)
        print(divHost)
        # get the pool size
        poolsize = len(divHost)
        print(poolsize)
        # dump the hostfiles
        self.hostInfo = self.dumpHost(divHost)
        return self.hostInfo, poolsize, totalproc

    def getProc(self):
        hostdict = {}
        if self.setmasternode != 0:
            # first = True
            print('set masternode')
        # else:
        #     first = False

        # get environment variable of host file list
        hostfile = os.environ["PE_HOSTFILE"]
        print(hostfile)
        fp = open(hostfile, "r")
        for x in fp:
            line = x.split()
            hostdict[line[0]] = int(line[1])
        fp.close()

        #        hostfile = os.environ["PBS_NODEFILE"]
        #        print hostfile
        #        fp = open(hostfile,"r")
        #        for x in fp:
        #            if first:
        #                masternode =x.split()[0]
        #                first= False
        #                continue
        #            line = x.split()
        #            hostdict[line[0]] = 12 #int(line[1]) zpliu
        #        fp.close()
        #        if self.setmasternode: hostdict.pop(masternode)
        return hostdict

    # chuck the list
    def chunks(self, arr, n):
        all = len(arr) - len(arr) % n
        return [tuple(arr[i:i + n]) for i in range(0, all, n)]

    def alloProc(self, hostdict, size=1):
        # construct a list of possible hosts
        totProc = 0
        hostList = []
        for key, val in hostdict.items():
            totProc += val
            hostList.extend([key] * val)

        print("Availiable proc number: ", totProc)
        if size == 0:
            size = totProc
        return self.chunks(hostList, size), totProc

    # dump hosts info into files
    # return a list of file names and corresponding host number

    def dumpHost(self, divHost):
        os.chdir(self.workdir)
        hostInfo = []
        ft = open(".hostfile", "w")
        for i, record in enumerate(divHost):
            fn = ".hostfile_%03d" % i
            fp = open(fn, "w")
            for line in record:
                fp.write(line + "\n")
                ft.write(line + "\n")
            fp.close()
            # why plus one ? in fact, I don't know either
            hostInfo.append((fn, len(set(record))))
        ft.close()
        return hostInfo


class allstr(list, Str):
    def __init__(self):
        list.__init__(self)

    def readfile(self, inputfile):
        f = open(inputfile, 'r')
        currentStr = -1
        for line in f:
            if ('Energy' in line or 'React' in line or 'Str' in line):
                self.append(Str())
                currentStr += 1
                try:
                    self[currentStr].energy = float(line.split()[-1])
                    if self[currentStr].energy.is_integer():
                        self[currentStr].energy = float(line.split()[-2])
                except Exception:
                    self[currentStr].energy = float(line.split()[-2])
            elif 'CORE' in line:
                self[currentStr].addatom(line, 1)
            elif ('PBC' in line) and ('ON' not in line):
                self[currentStr].abc = [float(x) for x in line.split()[1:]]
        f.close()
        self.nstr = currentStr + 1
        for str in self:
            str.sortatombyele()
            str.calAtomnum()
            str.abc2lat()
            # str.screenuppersurf()

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

    def printall(self, outfile, marker=False):
        f = open(outfile, 'w')
        f.write('!BIOSYM archive 2\n')
        f.write('PBC=ON\n')
        istr = 0

        for str in self:
            istr = istr + 1

            if not marker:
                marker = '%d' % istr
            f.write('     Energy       %6d   %6d   %14.8f\n' %
                    (istr, istr, str.energy))
            f.write('!DATE  %s\n' % marker)
            f.write('PBC %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n' %
                    (str.abc[0], str.abc[1], str.abc[2], str.abc[3],
                     str.abc[4], str.abc[5]))
            # f.write('PBC %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n'%(100,100,100,90,90,90))
            for i, atom in enumerate(str.atom):
                f.write(
                    '%-3s %15.9f %15.9f %15.9f CORE %4d %2s %2s %8.4f %4d\n' %
                    (atom.elesymbol, atom.xyz[0], atom.xyz[1], atom.xyz[2], i +
                     1, atom.elesymbol, atom.elesymbol, atom.charge, i + 1))
            f.write('end\nend\n')

    def printList(self, list, outfile):
        f = open(outfile, 'w')
        f.write('!BIOSYM archive 2\n')
        f.write('PBC=ON\n')
        istr = 0

        for istr in list:
            str = self[istr]

            f.write('     Energy       %6d   %6d   %14.8f\n' %
                    (istr, istr, str.energy))
            f.write('!DATE %8d\n' % istr)
            f.write('PBC %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n' %
                    (str.abc[0], str.abc[1], str.abc[2], str.abc[3],
                     str.abc[4], str.abc[5]))
            # f.write('PBC %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n'%(100,100,100,90,90,90))
            for i, atom in enumerate(str.atom):
                f.write(
                    '%-3s %15.9f %15.9f %15.9f CORE %4d %2s %2s %8.4f %4d\n' %
                    (atom.elesymbol, atom.xyz[0], atom.xyz[1], atom.xyz[2], i +
                     1, atom.elesymbol, atom.elesymbol, atom.charge, i + 1))
            f.write('end\nend\n')


class ReactionPair(allstr):
    def __init__(self, file):
        list.__init__(self)
        self.readfile(file)
        self.conn = 0
        self.TSstr = allstr()
        self.totdis = 0
        self.maxdis = 0
        return

    def buildfolder(self, dir):
        self.splituncm()
        os.system('mv path-0 %s' % dir)
        return

    def optfolder(self):
        maxP = os.popen('grep IS/M lasp.out').readline().strip()
        re77 = allstr()
        re77.readfile('SSWPath.arc')
        opt = allstr()
        opt.append(re77[0])
        maxP = max(re77, key=lambda x: x.energy)
        opt.append(maxP)
        opt.append(re77[-1])
        optfolder = opt.singlesplit()
        for dir in optfolder:
            os.system("sed -i 's/Run_type.*/Run_type 5/g' %s/input" % dir)
            os.system("sed -i 's/SSW.SSWsteps.*/SSW.SSWsteps 1/g' %s/input" %
                      dir)
        return optfolder

    def checkconn(self, dirs):
        result = []
        os.system('cat %s/all.arc > check.arc' % (dirs[0]))
        length = int(os.popen('cat check.arc | wc -l').readline().strip())
        os.system('tail -n%d %s/all.arc >> check.arc' % (length - 2, dirs[-1]))
        lsame, name1, name2 = checksame('check.arc')
        if (self[0].name != name1) or (self[1].name != name2):
            print('wrong pair')
            # return self
        for i in range(len(dirs) - 1):
            os.system('cat %s/all.arc > check.arc' % (dirs[i]))
            os.system('tail -n%d %s/all.arc >> check.arc' %
                      (length - 2, dirs[i + 1]))
            lsame, name1, name2 = checksame('check.arc')
            if lsame == 0:
                tmp = ReactionPair('check.arc')
                tmp[0].name = name1
                tmp[1].name = name2
                tmp.ISFSdistance()
                result.append(tmp)
        return result

    def extra(self):
        length = int(os.popen('cat TSmode.arc | wc -l').readline().strip())
        if length < 6:
            # print 'can not locate TS'
            return -1
        TSmode = allstr()
        TSmode.readfile('TSmode.arc')
        extra = allstr()
        extra.append(TSmode[0])
        extra.append(TSmode[-1])
        extra.splituncm()
        os.system('mv path-0 extra')
        return 1

    def checkextra(self):
        length = int(os.popen('cat SSWPath.arc | wc -l').readline().strip())
        if length < 6:
            # print 'can not connect in extra'
            return -1, self
        re77 = allstr()
        re77.readfile('SSWPath.arc')
        extra = allstr()
        extra.append(re77[0])
        extra.append(re77[-1])
        extra.printall('extra.arc')
        extra = ReactionPair('extra.arc')

        # simplly for len(self) == 2
        for i in range(len(self)):
            check = allstr()
            check.append(self[i])
            check.append(extra[i])
            check.printall('check.arc')

            result, self[i].name, extra[i].name = checksame('check.arc')

            check = allstr()
            check.append(self[i])
            check.append(extra[1 - i])
            check.printall('check.arc')
            result2, self[i].name, extra[1 - i].name = checksame('check.arc')

            if (result == 1 or result2 == 1):
                continue
            else:
                # print 'fail to extra'
                return 0, self
        try:
            extra.TS = float(
                os.popen('grep IS/TS/FS lasp.out').readline().split()[3])
        except Exception:
            return 0, self
        if (extra.TS - extra[0].energy) < 0 or (extra.TS -
                                                extra[1].energy) < 0:
            return 0, self
        extra.ISFSdistance()
        if (extra.totdis > (self.totdis + 0.3)
                or abs(extra.TS - self.TS) > 0.2):
            return 0, self
        TSstr = allstr()
        TSstr.readfile('TSstr.arc')
        extra.append(TSstr[0])
        extra.TSstr.readfile('TSstr.arc')
        return 1, extra

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


class Path(ReactionPair):
    def __init__(
        self,
        dir,
    ):
        list.__init__(self)
        self.dir = dir
        self.append(ReactionPair('%s/uncm.arc' % self.dir))

    def run(self, depth, runmode, *argv):
        # print argv
        os.chdir(self.dir)
        ipath = int((self.dir).split('-')[-1])
        fconn = open('conninfo', 'w')
        self[0].ISFSdistance()
        _tmp, self[0][0].name, self[0][1].name = checksame('uncm.arc')
        if ('94m' in self[0][0].name) or ('94m' in self[0][1].name):
            fconn.write('wrong pair\n')
            fconn.close()
            return self

        if (self[0][0].name == self[0][1].name):
            fconn.write('wrong pair: same IS/FS\n')
            fconn.close()
            return self

        fconn.write('initpair  %14.8f   %14.8f   %12.8f  %s %s   conn:%d\n' %
                    (self[0][0].energy, self[0][1].energy, self[0].totdis,
                     self[0][0].name, self[0][1].name, self[0].conn))

        for icycle in range(depth):

            fconn.write('cycle-%d\n' % icycle)
            os.chdir(self.dir)
            print('start cycle %d' % icycle)
            result = self.runcycle(icycle, runmode, *argv)

            for i, pair in enumerate(self):
                if ('94m' in pair[0].name) or ('94m' in pair[1].name):
                    pair.conn = -3

                fconn.write(
                    'sect-%d    %14.8f   %14.8f   %12.8f  %s %s   conn:%d\n' %
                    (i, pair[0].energy, pair[1].energy, pair.totdis,
                     pair[0].name, pair[1].name, pair.conn))
#            scount = 0
            if result != 0:
                break

#            for pair in self:
#                if pair.conn ==1:
#                    scount=scount +1
#            if scount == len(self):
#                result=1
#                break
        fconn.close()
        if result != 1:
            # print 'fail to connect string'
            os.chdir(self.dir)
            fout = open('pathinfo', 'w')
            fout.write('fail to connect string\n')
            fout.close()

        else:
            # print ('TS ')
            lTS = self.findTS(runmode, *argv)
            os.chdir(self.dir)
            # print ('dump final path')
            fout = open('pathinfo', 'w')
            fgoodsect = open('goodsect', 'w')
            for i, pair in enumerate(self):
                marker = '%s' % (self.dir) + '/TS/sect-%d/' % (i)
                pair.printall('sect-%d.arc' % i, marker)
                os.system('cat sect-%d.arc >> final.arc' % i)
                try:
                    pair.barrier = pair.TS - pair[0].energy
                    fout.write(
                        'path-%d:sect-%d IS/TS/FS    %14.8f %14.8f %14.8f %12.8f %12.8f  %s %s\n'
                        % (ipath, i, pair[0].energy, pair.TS, pair[1].energy,
                           pair.totdis, pair.barrier, pair[0].name,
                           pair[1].name))
                except Exception:
                    pair.barrier = pair.maxP - pair[0].energy
                    fout.write(
                        'path-%d:sect-%d IS/maxP/FS  %14.8f %14.8f %14.8f %12.8f %12.8f  %s %s\n'
                        % (ipath, i, pair[0].energy, pair.maxP, pair[1].energy,
                           pair.totdis, pair.barrier, pair[0].name,
                           pair[1].name))
                if pair.lTS == 1:
                    os.system('cat sect-%d.arc >> goodsect.arc' % i)
                    pair.barrier = pair.TS - pair[0].energy
                    fgoodsect.write(
                        'path-%d:sect-%d IS/TS/FS    %14.8f %14.8f %14.8f %12.8f %12.8f  %s %s\n'
                        % (ipath, i, pair[0].energy, pair.TS, pair[1].energy,
                           pair.totdis, pair.barrier, pair[0].name,
                           pair[1].name))

            fout.close()
            fgoodsect.close()

            if lTS == 0:
                os.system('cp pathinfo goodpath')
        return self

    def runcycle(self, icycle, runmode, *argv):
        os.mkdir('cycle-%d' % icycle)
        cycledir = os.path.join(self.dir, 'cycle-%d' % icycle)
        os.system("sed -i 's/Run_type.*/Run_type 2/g'  input")
        os.system("sed -i 's/DESW.task.*/DESW.task string/g' input")

        os.system('cp input cycle-%d' % icycle)
        os.system('ln -s  %s/*.pot cycle-%d' % (os.getcwd(), icycle))
        result = {}
        # print self[0][0].energy
        for i, pair in enumerate(self):
            os.chdir(cycledir)
            pair.printall('sect-%d.arc' % i)
            if (pair.conn != 0):
                continue

            print('start cal sect-%d' % i)
            pair.buildfolder('sect-%d' % i)
            workdir = os.path.join(cycledir, 'sect-%d' % i)
            if runmode == 'local':
                runprog_local(workdir, *argv)
            elif runmode == 'cluster':
                runprog_cluster(workdir, *argv)

            length = int(
                os.popen('cat SSWPath.arc | wc -l').readline().strip())
            if length < 8:
                pair.conn = -1
                return -1
            # pair.totdis= float(os.popen('grep distance lasp.out').readline().split()[3])
            pair.maxP = float(
                os.popen('grep IS/MaxP/FS lasp.out').readline().split()[3])
            print(pair.totdis)
            optfolder = pair.optfolder()
            for dir in optfolder:
                optdir = os.path.join(workdir, dir)
                if runmode == 'local':
                    runprog_local(optdir, *argv)
                elif runmode == 'cluster':
                    runprog_cluster(optdir, *argv)
            os.chdir(workdir)
            tmpresult = pair.checkconn(optfolder)
            if len(tmpresult) == 1:
                # tmpresult[0].ISFSdistance()
                if tmpresult[0].totdis < (pair.totdis - 0.3):
                    result[i] = tmpresult
                else:
                    pair.conn = 1
                    # if not hasattr(pair[0], 'name'):
                    #    pair.printall('check.arc')
                    #    tmpre,pair[0].name,pair[1].name = checksame('check.arc')
                    result[i] = [pair]
            else:
                for newpair in tmpresult:
                    if newpair.totdis > (pair.totdis + 0.3):
                        pair.conn = -2
                        return -2
                result[i] = tmpresult
        pcount = 0
        for key, val in result.items():
            print(val)
            if len(val) == 1:
                abspos = key + pcount
                self.pop(abspos)
                self.insert(abspos, val[0])
            else:
                abspos = key + pcount
                self.pop(abspos)
                for i, pair in enumerate(val):
                    self.insert(i + abspos, pair)
                pcount = pcount + len(val) - 1

#        if pcount == 0:
#            return 1
#        else :
#            return 0

        for pair in self:
            if pair.conn != 1:
                return 0
        return 1

    def findTS(self, runmode, *argv):
        os.chdir(self.dir)
        finfo = open('TSinfo', 'w')
        os.mkdir('TS')
        pwd = os.path.join(self.dir, 'TS')
        os.chdir(pwd)
        os.system('cp ../input .')
        os.system('ln -s  ../*.pot .')
        os.system("sed -i 's/Run_type.*/Run_type 2/g'  input")
        os.system("sed -i 's/DESW.task.*/DESW.task TS/g' input")

        for i, pair in enumerate(self):
            os.chdir(pwd)
            pair.buildfolder('sect-%d' % i)
            workdir = os.path.join(pwd, 'sect-%d' % i)
            if runmode == 'local':
                runprog_local(workdir, *argv)
            elif runmode == 'cluster':
                runprog_cluster(workdir, *argv)
            lTS = pair.extra()
            if lTS != 1:
                finfo.write('sect-%d: can not locate TS\n' % i)
                pair.lTS = -2
                continue

            pair.TS = float(
                os.popen('grep IS/TS/FS lasp.out').readline().split()[3])

            workdir = os.path.join(workdir, 'extra')
            if runmode == 'local':
                runprog_local(workdir, *argv)
            elif runmode == 'cluster':
                runprog_cluster(workdir, *argv)
            # print 'checkextra'
            result, extrapair = pair.checkextra()
            # print result,len(extrapair)
            pair.lTS = result
            if result == 1:
                finfo.write('sect-%d: TS search success\n' % i)
                # extrapair.ISFSdistance()
                extrapair.lTS = result
                self.pop(i)
                self.insert(i, extrapair)
            elif result == -1:
                finfo.write('sect-%d: can not connect in extra\n' % i)
            elif result == 0:
                finfo.write('sect-%d: fail to extra\n' % i)

            # TSstr=allstr()
            # TSstr.readfile('../TSstr.arc')
            # pair.append(TSstr[0])
            # air.TSstr.readfile('../TSstr.arc')

        for pair in self:
            if pair.lTS != 1:
                finfo.close()
                return -1
        finfo.write('all success\n')
        finfo.close()
        return 0


#    def optfolder(self):
#        maxP =os.popen('grep IS/M output').readline().strip()
#        re77=allstr()
#        re77.readfile('re77.arc')
#        opt = allstr()
#        opt.append(re77[0])
#        maxP =max(re77,key=lambda x: x.energy)
#        opt.append(maxP)
#        opt.append(re77[-1])
#        optfolder= opt.singlesplit()
#        for dir in optfolder:
#            os.system("sed -i 's/Run_type.*/Run_type 5/g' %s/input"%dir)
#            os.system("sed -i 's/SSW.SSWsteps.*/SSW.SSWsteps 1/g' %s/input"%dir)
#        return optfolder

#    def checkconn(self,dirs):
#        result= []
#        for i in range(len(dirs)-1):
#            os.system('cat %s/all.arc > check.arc'%(dirs[i]))
#            length = int(os.popen('cat check.arc | wc -l').readline().strip())
#            os.system('tail -n%d %s/all.arc >> check.arc'%(length-2,dirs[i+1]))
#            lsame, name1,name2= checksame('check.arc')
#            if lsame == 0:
#                tmp = ReactionPair('check.arc')
#                tmp[0].name = name1
#                tmp[1].name = name2
#                tmp.ISFSdistance()
#                result.append(tmp)
#        return result
#
#
#    def extra(self):
#        TSmode=allstr()
#        TSmode.readfile('TSmode.arc')
#
#        return


class allPath(Path):
    def __init__(self, dir, prog, workdirs, depth):
        list.__init__(self)
        self.dir = dir
        self.prog = prog
        self.depth = depth
        self.workdirs = workdirs
        for workdir in workdirs:
            absdir = os.path.join(self.dir, workdir)
            self.append(Path(absdir))

    def runall_local(self, poolsize, ncpu, runmode):
        # Ncore= 72
        # poolsize= Ncore/ncpu
        pool = Pool(processes=poolsize)
        # result = []
        for path in self:
            print(path.dir)
            pool.apply_async(run_wrapper,
                             args=(path, self.depth, runmode, self.prog, ncpu))
            # path.run(depth,self.prog,ncpu)
        pool.close()
        pool.join()

        self.collectinfo()
        return

    def runall_cluster(self, poolsize, ncpu, hostinfo, runmode, poolcount=0):
        pool = Pool(processes=poolsize)
        # result = []
        for path in self:
            print(path.dir)
            pool.apply_async(run_wrapper,
                             args=(path, self.depth, runmode, self.prog, ncpu,
                                   hostinfo, self.dir, os.environ, poolcount))
        pool.close()
        pool.join()

        self.collectinfo()
        return

    def collectinfo(self):
        os.chdir(self.dir)
        os.system('rm -f allpath')
        os.system('rm -f allgoodpath')
        os.system('rm -f allgoodsect')
        os.system('rm -f allgoodsect.arc')
        for workdir in self.workdirs:
            os.system('echo %s >> allpath' % workdir)
            os.system('cat %s/pathinfo >> allpath' % workdir)
            os.system('cat %s/goodpath >> allgoodpath' % workdir)
            os.system('cat %s/goodsect >> allgoodsect' % workdir)
            os.system('cat %s/goodsect.arc >> allgoodsect.arc' % workdir)


def run_wrapper(instance, *argv):
    return instance.run(*argv)


def checksame(file):
    #    name2= result[1].split()[0]
    result = os.popen('minimum_new -i %s' % file).readlines()
    name1 = reduce(lambda a, b: a + b,
                   ["%s" % s for s in result[3].split()[3:]])
    name2 = reduce(lambda a, b: a + b,
                   ["%s" % s for s in result[4].split()[3:]])
    #    if name1 not in namedict.keys():
    #
    #        self.namedict = {}
    #        for i in range(len(self)):
    #            if self[i].upper ==1: continue
    #            name = self[i].name
    #            if name not in self.namedict.keys():
    #               self.namedict[name]=[i,]
    #            else:
    #               self.namedict[name].append(i)
    #        self.nstrpername =[(x,len(y)) for x,y in self.namedict.items()]
    #        self.nstrpername.sort(key =lambda X:X[1], reverse = True)

    if name1 == name2:
        lsame = 1
    else:
        lsame = 0
    return lsame, name1, name2


def checksame2(file):
    check = allstr()
    check.readfile(file)
    lsame = 0
    if abs(check[2].energy - check[2].energy) < 0.05:
        lsame = 1
    return lsame


def runprog_local(workdir, prog, ncpus):
    try:
        os.chdir(workdir)
        mpiprog = "mpirun -np %d " % (ncpus) + prog
        fout = open('output', 'w')
        subprocess.call(mpiprog,
                        stdout=fout,
                        stderr=fout,
                        shell=True,
                        executable='/bin/bash')
        fout.close()
        return
    except Exception as e:
        traceback.print_exc()
        raise e
    # except:
    #     print('error')
    #     return


def runprog_cluster(workdir, prog, ncpus, hostInfo, rootdir, env, poolcount=0):

    try:
        os.chdir(rootdir)
        cwd = os.getcwd()
        # exit = False
        print(workdir)
        print(multiprocessing.current_process().name)
        nodeInfo = hostInfo[
            int(multiprocessing.current_process().name.split("-")[-1]) - 1 -
            poolcount]
        os.chdir(workdir)
        mf = os.path.join(cwd, nodeInfo[0])
        # mpiprog = "/home/software/mpi/intel/impi/4.0.1.007/bin64/mpirun --rsh=ssh -machinefile %s -np %d "%(mf, ncpus) + prog
        mpiprog = "mpirun -r ssh -machinefile %s -np %d " % (mf, ncpus) + prog
        fout1 = open("proginfo", "w")
        fout1.write("Current process: " +
                    multiprocessing.current_process().name + "\n")
        fout1.write("Host file: " + mf + "\n")
        fout = open('output', 'w')
        totalrun = 'source ~/.bashrc; ' + mpiprog
        # nodename = os.popen('head -1 %s'%mf).readline().strip()
        # totalrun = 'source ~/.bashrc; ssh '+nodename +"; "+mpiprog
        fout1.write(totalrun + "\n")
        fout1.close()
        # child= subprocess.Popen(mpiprog,stdout = fout, stderr=fout,shell=True,executable='/bin/bash',preexec_fn = os.setpgrp)
        subprocess.call(mpiprog,
                        stdout=fout,
                        stderr=fout,
                        shell=True,
                        executable='/bin/bash')
        # print ('start run job in %s'%workdir)
        # fout1.write('pid   %d\n'%child.pid)
        # fout1.close()
        #        pid =child.pid
        #        if maxtime:
        #            alltime = 0
        #        while not exit:
        #            time.sleep(3)
        #            returnCode= child.poll()
        #            if glob.glob('killsignal'):
        #                os.kill(-pid,9)
        #                time.sleep(3)
        #                fout.write('kill %s\n'%pid)
        #                a=os.waitpid(pid,0)
        #                print a
        #                time.sleep(3)
        #                exit =True
        #            if isinstance(returnCode,int):
        #                if returnCode == 0:
        #                    fout.write('successfully done\n')
        #                else:
        #                    fout.write('something wrong: returnCode  %d\n'%returnCode)
        #                exit =True
        #            if maxtime:
        #                alltime = alltime+30
        #                if alltime >maxtime:
        #                    os.kill(-pid,9)
        #                    time.sleep(3)
        #                    fout.write('time out :kill %s\n'%pid)
        #                    a=os.waitpid(pid,0)
        #                    print a
        #                    time.sleep(3)
        #                    exit =True
        #
        #        #fout1.write(str(child.poll())+'\n')
        #        fout.write('exit\n')
        fout.close()
        return
    except Exception as e:
        traceback.print_exc()
        raise e
    # except:
    #     print('error')
    #     return


if __name__ == "__main__":

    #    fout = open('output','w')
    #    subprocess.call('mpirun -np 12 /home10/bin/lasp-release-1.0/lasp',stdout = fout, stderr=fout,shell=True,executable='/bin/bash')
    #    fout.close()
    #
    os.system("sed -i 's/SSW.Pathway.*/SSW.Pathway F/g' input")
    rootdir = os.getcwd()

    os.system('rm -rf path-*')

    alluncm = allstr()
    alluncm.readfile('uncm.arc')
    workdirs = alluncm.splituncm()
    i = max(alluncm, key=lambda x: x.energy)

    print(sys.argv)

    prog = '/home10/bin/lasp-release-1.0/lasp-beta5'
    #    prog ='/home10/bin/sswoop-main.current'
    ncpu = 4
    depth = 4
    allpath = allPath(rootdir, prog, workdirs, depth)

    # global namedict
    # namedict={}

    runmode = sys.argv[1]
    if runmode == 'local':
        print('local mode')

        Ncore = 72
        poolsize = Ncore / ncpu
        allpath.runall_local(poolsize, ncpu, runmode)
    elif runmode == 'cluster':
        print('cluster mode')

        Host = Hostfile(rootdir, ncpu, 0)
        hostInfo, poolsize, totalproc = Host.setHostfile()
        allpath.runall_cluster(poolsize, ncpu, hostInfo, runmode)
