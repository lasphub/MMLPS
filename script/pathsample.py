#! /usr/bin/env /home10/kpl/anaconda2/envs/rdkitenv/bin/python

# from rdkit import Chem
# from rdkit.Chem import AllChem
# from rdkit import DataStructs
import re
import os
import glob
from multiprocessing import Pool
import time
import subprocess
import numpy as np

from allstr_uncm import AllStr
from restart import runprog_local
from hostfile_tmp import Hostfile, runprog_cluster
from analyze_surfacemode import AllPair
from allpair import AllPair as AllPair_Screen
from ECFP import AllStr as AllStr_ECFP
from pathfindingYen import PathanalyzeYen
from endCheckMK import EndCheckMK
from pathsample_para import pathSampleMultiPara, pathSamplePara

direcToBallance = pathSampleMultiPara["direcToBallance"]
direcJoint = pathSampleMultiPara["direcJoint"]


class DElink(AllStr):
    def __init__(self, file, rootdir, pathSamplePara):
        AllStr.__init__(self)
        self.readfile(file)
        # self.GetAllName(file)
        self.rootdir = rootdir
        self.Lrestart = pathSamplePara["restart"]
        self.ncpusample = pathSamplePara["ncpusample"]
        self.ncpusplit = pathSamplePara["ncpusplit"]
        self.MaxNumofSinglePair = pathSamplePara["MaxNumofSinglePair"]
        self.ncputotal = pathSamplePara["ncputotal"]
        self.Lscreenupper = pathSamplePara["Lscreenupper"]
        self.ncputotal = pathSamplePara["ncputotal"]

        _tmp = AllStr_ECFP()
        _tmp.readfile(file)
        # _tmp.GetAllpuresminame()
        # _tmp.GenCanoicalSmiles()
        # self.namelist =[Str.canname for Str in _tmp]

        # surfacemode
        _tmp.GetAllSminame_FromEcfp(numproc=24, colorflag=0)
        _tmp.GetAllECFPname(numproc=24)
        self.namelist = [Str.ECFPname for Str in _tmp]
        for struc1, struc2 in zip(self, _tmp):
            struc1.ECFPname = struc2.ECFPname
            struc1.sminame = struc2.sminame

        self.prog = pathSamplePara["prog"]
        self.runmode = pathSamplePara["runmode"]
        self.poolcount = 96

        self.pairdata = AllPair()

        self.printlist([0], 'target.arc')
        if not glob.glob('analyze'):
            os.mkdir('analyze')
        self.analyzedir = os.path.join(self.rootdir, 'analyze')
        os.system('cp target.arc analyze')
        self.child = subprocess.Popen('echo test > subprocesstest',
                                      shell=True,
                                      executable='/bin/bash',
                                      preexec_fn=os.setpgrp)

        for Str in self:
            Str.sample = self.Lrestart
        self.updateinit = 1
        self.Gendatamode = 0

    def Run(self):

        print("start sampling ...")
        if (self.Lrestart == 1):
            self.SelectStr_new()

        for icycle in range(0, 100):
            os.chdir(self.rootdir)
            if glob.glob('cycle_%d' % icycle):
                continue
            self.printall('all_cycle_%d.arc' % icycle)
            endFlag = self.RunCycle(icycle)
            if endFlag is True:
                print("End!")
                os.system(
                    """if [ "$USER" = "shiyf" ]; then iyuuMsg.py -jd; fi""")
                break
        return

    def UppdateInit(self):
        tmp = AllStr_ECFP()
        tmp.readfile('para0/best.arc')
        # tmp.GetAllsminame()
        tmp.GetTmpFakebmx()
        tmp.ScreenUpper()
        self.poolcount = self.poolcount + 48
        if len(tmp) == 0:
            return
        elif tmp[-1].upper == 1:
            return
        else:
            tmp[-1].GenECFPname()
            if self.namelist[0] != tmp[-1].ECFPname or self.updateinit == 1:
                self.pop(0)
                self.insert(0, tmp[-1])
                self[0].sample = 1
                self.namelist[0] = tmp[-1].ECFPname
                tmp.printlist([-1], 'target.arc')
                os.system('cp target.arc ../analyze/')
            return

    def RunCycle(self, icycle):
        os.chdir(self.rootdir)
        # if glob.glob('cycle_%d'%icycle):
        #    return
        os.mkdir('cycle_%d' % icycle)
        os.chdir('cycle_%d' % icycle)

        # rootdir = os.path.join(self.rootdir, 'cycle_%d' % icycle)
        os.system('cp ../sourcedir/*.pot .')
        os.system('cp ../sourcedir/input .')
        os.system('cp ../sourcedir/input_split .')

        t0 = time.time()
        self.PathSample()

        if len(self) == 1:
            self.UppdateInit()

        t1 = time.time()
        print('Sample time: %f' % (t1 - t0))

        self.CalPath_newprog()
        t2 = time.time()
        print('calpath time: %f' % (t2 - t1))

        # SelectStr new

        nadd = self.SelectStr_new()
        t3 = time.time()
        print('select time: %f' % (t3 - t2))

        endFlag = self.EndCheck()
        t4 = time.time()
        print('EndCheck time: %f' % (t4 - t2))

        if nadd == 0 or endFlag is True:
            return True
        else:
            return False

    def collectalllow(self):
        os.system('cat para*/best.arc >> allbest.arc')
        tmp = AllStr_ECFP()
        tmp.readfile('allbest.arc')

        if len(tmp) > 0:
            ibest = -1
            tmp.ScreenUpper()
            tmp.sort(key=lambda X: X.energy)
            for iprint in range(len(tmp)):
                if tmp[iprint].upper == 0:
                    ibest = iprint
                    break
            if ibest != -1:
                tmp.printlist([ibest], 'best.arc')
                os.system('cat best.arc >>  ../analyze/addmin.arc')
            else:
                os.system('cat sample.arc >>  ../analyze/addmin.arc')

        else:
            os.system('cat sample.arc >>  ../analyze/addmin.arc')

    def PreSelectPair(self, nmax=1000, nmaxsingle=150, Lscreenupper=0):
        _tmp = AllPair_Screen()
        _tmp.readfile("uncm.arc", Lscreenupper)
        _tmp.OutPair(nmax, nmaxsingle, "outpair.arc")

    def GetSampledStruc(self):

        sampledStruc = AllStr_ECFP()
        for direc, ballanceName in direcToBallance.items():
            allCycleArc = glob.glob(direc + '/all_cycle*.arc')
            if len(allCycleArc) == 0:
                continue
            allCycleArc = sorted(
                allCycleArc, key=lambda x: int(re.search(r"\d+", x).group()))
            _tmp = AllStr_ECFP()
            _tmp.readfile(allCycleArc[-1])
            _tmp.GetAllSminame_FromEcfp(numproc=24, colorflag=0)
            _tmp.GetAllECFPname(numproc=24)
            self.poolcount += 96

            for struc in _tmp:
                fragments = [
                    re.sub(r"-.*ads", "", frag)
                    for frag in struc.sminame.split('.')
                ]
                _normSmiName = fragments + ballanceName.split(".")
                struc.normSminame = ".".join(sorted(_normSmiName))
            sampledStruc.extend(_tmp)
            print(direc)
            print([struc.normSminame for struc in _tmp])
        sampledNormSminame = [struc.normSminame for struc in sampledStruc]
        return sampledNormSminame

    def SelectStr_new(self):
        """for first cycle"""
        os.chdir(self.analyzedir)

        os.system('cat allgoodsect.arc_tmp >> allgoodsect.arc')
        os.system('cat allgoodsect_tmp >> allgoodsect')
        os.system('rm -f allgoodsect.arc_tmp')
        os.system('rm -f allgoodsect_tmp')

        sampledNormSminame = self.GetSampledStruc()
        ballanceName = direcToBallance.get(self.rootdir, "")
        name, nextStr = PathanalyzeYen(self.namelist[0], sampledNormSminame,
                                       ballanceName)
        if name is None:
            return 0
        self.poolcount += 96 * 2 + 24
        nextStr.sample = 0
        self.append(nextStr)
        self.namelist.append(name)
        os.chdir(self.rootdir)
        return 1

    def PathSample(self):
        rootdir = os.getcwd()

        _cyclestr = AllStr()
        for Str in self:
            if Str.sample == 0:
                _cyclestr.append(Str)
                Str.sample = 1
        _cyclestr.printall('sample.arc')
        _cyclestr.singlesplit()

        if self.runmode == 'local':
            Ncore = self.ncputotal
            ncpu = self.ncpusample
            poolsize = int(Ncore / ncpu)

            # if poolsize > len(_cyclestr):
            #    os.system('cp ')

            pool = Pool(processes=poolsize)
            result = []
            njob = 1 * poolsize
            for i in range(njob):
                if i > 0:
                    os.system('cp -r para0 para%d' % i)
                path = os.path.join(rootdir, 'para%d' % i)
                result = pool.apply_async(runprog_local,
                                          args=(path, self.prog, ncpu))
            pool.close()
            pool.join()

        elif self.runmode == 'cluster':
            ncpu = self.ncpusample
            Host = Hostfile(rootdir, ncpu, setmasternode=0)
            if glob.glob('../hostfile_now'):
                hostInfo, poolsize, totalproc = Host.setHostfile(
                    '../hostfile_now')
            else:
                hostInfo, poolsize, totalproc = Host.setHostfile()

            pool = Pool(processes=poolsize)
            result = []
            njob = poolsize

            for i in range(njob):
                if i > 0:
                    os.system('cp -r para0 para%d' % i)
                path = os.path.join(rootdir, 'para%d' % i)
                result = pool.apply_async(runprog_cluster,
                                          args=(path, self.prog, ncpu,
                                                hostInfo, rootdir,
                                                self.poolcount))
            # workdir, prog, ncpus, hostInfo, rootdir, env, poolcount=0
            pool.close()
            pool.join()
            self.poolcount = self.poolcount + poolsize

        os.chdir(rootdir)
        for i in range(njob):
            os.system('cat para%d/uncm.arc >> uncm.arc' % (i))

        self.collectalllow()

    def CalPath_newprog(self):
        self.PreSelectPair(self.ncputotal, self.MaxNumofSinglePair,
                           self.Lscreenupper)
        pwd = os.getcwd()

        allpair = AllStr()
        allpair.readfile('outpair.arc')
        split = int(len(allpair) / 2000) + 1
        for i in range(split):
            if i != split - 1:
                allpair.printlist(range(2000 * (i), 2000 * (i + 1)),
                                  'uncm.arc_%d' % i)
            else:
                allpair.printlist(range(2000 * (i), len(allpair)),
                                  'uncm.arc_%d' % i)

        for isplit in range(split):
            os.chdir(pwd)
            if glob.glob('split_%d' % isplit):
                continue
            os.mkdir('split_%d' % isplit)
            os.chdir('split_%d' % isplit)
            rootdir = os.getcwd()

            os.system('cp ../uncm.arc_%d .' % isplit)
            os.system('cp ../input_split input')
            os.system('cp ../*.pot .')

            allpair = AllStr()
            allpair.readfile('uncm.arc_%d' % isplit)
            workdirs = allpair.splituncm()
            print(workdirs)

            if self.runmode == 'local':
                Ncore = self.ncputotal
                ncpu = self.ncpusplit
                poolsize = int(Ncore / ncpu)

                # if poolsize > len(_cyclestr):
                #    os.system('cp ')

                pool = Pool(processes=poolsize)
                result = []

                for workdir in workdirs:
                    path = os.path.join(rootdir, workdir)
                    result = pool.apply_async(runprog_local,
                                              args=(path, self.prog, ncpu))

                pool.close()
                pool.join()

            elif self.runmode == 'cluster':

                ncpu = self.ncpusplit
                Host = Hostfile(rootdir, ncpu, setmasternode=0)
                if glob.glob('../../hostfile_now'):
                    hostInfo, poolsize, totalproc = Host.setHostfile(
                        '../../hostfile_now')
                else:
                    hostInfo, poolsize, totalproc = Host.setHostfile()

                pool = Pool(processes=poolsize)
                result = []

                for workdir in workdirs:
                    path = os.path.join(rootdir, workdir)
                    result = pool.apply_async(runprog_cluster,
                                              args=(path, self.prog, ncpu,
                                                    hostInfo, rootdir,
                                                    self.poolcount))
                    # workdir, prog, ncpus, hostInfo, rootdir, env, poolcount=0
                    # result= pool.apply_async(runprog_cluster,args= (path,self.prog,ncpu,hostInfo,rootdir,os.environ,self.poolcount))
                pool.close()
                pool.join()
                self.poolcount = self.poolcount + poolsize

            os.chdir(rootdir)
            for path in workdirs:
                os.system('cat %s/RPTS.arc >> allgoodsect.arc' % path)
            os.system(
                'cat allgoodsect.arc >> ../../analyze/allgoodsect.arc_tmp')

    def EndCheck(self):
        if not os.path.isdir(direcJoint):
            self.MakeJointDirec()
        if not hasattr(self, "oldProdRate"):
            self.oldProdRate = 1e90
            self.endCountFlag = 0

        prodRate = EndCheckMK(direcToBallance.keys(),
                              [("CO2.H2.H2.H2", "CH3OH.H2O")],
                              direcJoint)
        if not prodRate:
            self.endCountFlag = 0
            self.oldProdRate = 1e90
        else:
            if abs(np.log(prodRate / self.oldProdRate)) < 1e-3:
                self.endCountFlag += 1
            else:
                self.endCountFlag = 0
            self.oldProdRate = prodRate

        if self.endCountFlag >= 8:
            endFlag = True
        else:
            endFlag = False
        return endFlag

    def MakeJointDirec(self):
        os.system("mkdir -p %s" % direcJoint)
        os.chdir(direcJoint)
        for direc in direcToBallance:
            direcName = direc.split("/")[-1]
            os.system("ln -s %s/analyze %s" % (direc, direcName))
        os.system('echo "startmole CO2_2/3H2\nstartsurf Cu111\n" > input')
        os.chdir(self.rootdir)


# def ReadPara(file):
#     f = open(file, 'r')
#     line = f.readline()
#     para = {}
#     while line:
#         if len(line.split()) == 0:
#             line = f.readline()
#             continue
#         elif line.split()[0] == 'totalthick':
#             para['totalthick'] = float(line.split()[1])
#         elif line.split()[0] == 'ncputotal':
#             para['ncputotal'] = int(line.split()[1])
#         elif line.split()[0] == 'ncpusample':
#             para['ncpusample'] = int(line.split()[1])
#         elif line.split()[0] == 'ncpusplit':
#             para['ncpusplit'] = int(line.split()[1])
#         elif line.split()[0] == 'Potlib':
#             para['Potlib'] = line.split()[1]
#         elif line.split()[0] == 'Bulklib':
#             para['Bulklib'] = line.split()[1]
#         elif line.split()[0] == 'Mollib':
#             para['Mollib'] = line.split()[1]
#         elif line.split()[0] == 'prog':
#             para['prog'] = line.split()[1]
#         elif line.split()[0] == 'restart':
#             para['restart'] = int(line.split()[1])
#         elif line.split()[0] == 'Lscreenupper':
#             para['Lscreenupper'] = int(line.split()[1])
#         elif line.split()[0] == 'job':
#             para['job'] = line.split()[1]
#         elif line.split()[0] == 'runmode':
#             para['runmode'] = line.split()[1]
#         elif line.split()[0] == 'MaxNumofSinglePair':
#             para['MaxNumofSinglePair'] = int(line.split()[1])
#         elif line.split()[0] == 'ThresholdofSelect':
#             para['ThresholdofSelect'] = int(line.split()[1])
#         # elif line.split()[0] == '%block':
#         #     blockname = line.split()[-1]
#         #     lines = []
#         #     blockline = f.readline().strip('\n')
#         #     while blockline.split()[0] != '%endblock':
#         #         lines.append(blockline)
#         #         blockline = f.readline().strip('\n')
#         # if blockname == 'ReactPattern':
#         #     para[blockname] = PatternBlock(lines)
#         # else:
#         #     para[blockname] = lines
#         line = f.readline()
#     f.close()
#     return para

if __name__ == "__main__":

    rootdir = os.getcwd()

    # para = ReadPara('console')
    # test = DElink('start.arc', rootdir, para['runmode'], para['prog'],
    #               para['restart'], para['ncpusample'], para['ncpusplit'],
    #               para['MaxNumofSinglePair'], para['ThresholdofSelect'],
    #               para['Lscreenupper'], para['ncputotal'])
    test = DElink('start.arc', rootdir, pathSamplePara)
    test.Run()
    #test.SelectStr_new()
    #test.EndCheck()
