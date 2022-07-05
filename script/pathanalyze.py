#!/home10/shiyf/bin/miniconda3/bin/python
"""
    examine pathway of MaxP/TS type
"""
# from systemop import getsysargs
import sys
# from Font import endMark, outGreen, outBlue  # , outPink, outYellow,
# from analyze_organic import AllPair
from analyze_surfacemode import AllPair as AllPair_surfacemode
import analyze_surfacemode
from allstr_new import allstr
# from reactpair import Pair as rppair
import time
from datetime import datetime
import os
import re
from dbclient import LaspDBPathSmpl
import globalVar
# import hNNCalPhonon

"""
filemode 读文件格式 IS/FS/TS 开1 IS/TS/FS 开2  default 1  读allgoodsect.arc
readallmin  额外加GM，开2，文件名addmin.arc  开1 allname.arc
surface 一定要写
pts 是否打印过渡态
-all  打印所有GM
-allpair 打印所有反应
-targetpair  xx.arc    找一对反应的所有例子
-p 和一个反应物相关的所有反应对 快的 约等于allpair  1个结构
-l xx.arc  连接一对反应  给xx.arc  两个结构
-s 和一个反应物相关的所有产物  很慢很占内存
-multi 多个analyze一同进行分析，目前计划先写link
-LConly 仅识别碳原子上的成键变化
-LGibbs 加吉布斯校正，需要改源码
"""


def getsysargs(para, outpa={}):
    readflag = False
    for i, item in enumerate(para):
        if item[0] == '-':  # parameters
            readflag = True
            key = item[1:]
            outpa[key] = []
            continue
        if readflag:
            outpa[key].append(item)
    return outpa


def printall(self, fname):
    if len(self) == 0:
        return

    if len(self[0].atom) != 0:
        oldprintall(self, fname)
        return

    with open(fname, "w") as fp:
        fp.write("!BIOSYM archive 2\nPBC=ON\n")
        for struc in self:
            fp.write(dbLoader.pathSmplDownLoadArc(struc))


def AnalyzerCalOrDown(dbLoader, analyzer, direc):
    if dbLoader.PathSmplUploadCheckExist(direc):
        print("EXISTED in Database! Downloading .....")
        dbLoader.PathSmplDownLoadInfo(analyzer, direc)
        analyzer.LGibbs = LGibbs
        analyzer.GibbsCorrectAll()
        analyzer.GenAllStrMin()
        print("Downloading Done.")
        analyzer.GenPairdict_byname()

        # if 211
        # _tmp = allstr()
        # for pair in analyzer:
        #     atChanged = [at for at in pair.TSstr.atom if at.xyz[2] > 10 and 1.8 < at.xyz[0] < 2.7]
        #     if len(atChanged) > 0:
        #         print(atChanged[0].xyz)
        #         _tmp.append(atChanged)
        # analyzer = [pair for pair in analyzer if pair not in _tmp]
        # _tmp.printall("popout.arc")
        # print("Pop out Unwanted Structures")

    else:
        print("NOT EXISTED in Database! Calculate name and etc...")
        analyzer.readfile(filename=filename,
                          filemode=filemode,
                          Lallmininput=Lallmin,
                          LGibbs=LGibbs)
        print("Name Calculation done.")
        analyzer.GenPairdict_byname()

        if 'zpe' in para.keys():
            print("Calculating ZPE of (lowest TS of pair and all LM)")
            hNNCalPhonon.CalPathSmplZpe(analyzer, numproc=16)
            print("Calculating ZPE Done.")
        analyzer.GibbsCorrectAll()
        print("Gibbs correction done")
        print("Uploading to Database ......")
        dbLoader.PathSmplUpload(analyzer, globalVar.globalDict["startsurf"],
                                globalVar.globalDict["startmole"])
        print("Uploading Done.")


para = getsysargs(sys.argv, outpa={})

if 'h' in para.keys() or '-help' in para.keys():
    print(
        '********************************************************************')
    print(
        '*%s                    Analysis of pair information                    %s*'
        % (outGreen, endMark))
    print(
        '********************************************************************')
    print(
        '  -p strid/arcfile  :%s output all possible pair which contains input str%s'
        % (outBlue, endMark))
    print(
        '  -s strid/arcfile  :%s output all possible intermediates connecting to input str (sorted by barrier)%s'
        % (outBlue, endMark))
    print(
        '  -l linkpath       :%s filename of target pair arc, or id of target pair%s'
        % (outBlue, endMark))
    print('  -all              :%s output all str contains in dataset %s' %
          (outBlue, endMark))
    print('  -allpair          :%s output all pairs contains in dataset %s' %
          (outBlue, endMark))
    print(
        '********************************************************************')
    exit()

filemode = 1
Lallmin = 0
LGibbs = -1
LConly = 0
filename = "allgoodsect.arc"
# syf easy operation
# if "descore" in para.keys():
#     with open("lowestpath") as fp:
#         lines = fp.readlines()
#     path = []
#     for li in lines:
#         if "lowest" in li and '-' in li:
#             path.append([
#                 re.split("-+", li)[2],
#                 str(int(float(re.split("-+", li)[6])))
#             ])
#     for p in path:
#         score = p.pop()
#         i = 0
#         while len(score) > 3:
#             if float(score[0:2]) == 0:
#                 break
#             p.append(int(score[0:3]) / 100.0)
#             score = score[3:]
#             i += 1
#         print("-".join(["%.2f" % float(s) for s in p]))
#     exit()
if 'LConly' in para.keys():
    LConly = int(para['LConly'][0])

if "db" in para.keys():
    oldprintall = analyze_surfacemode.allstr.printall
    analyze_surfacemode.allstr.printall = printall

# t1 = time.clock()
t1 = time.perf_counter()

if 'filemode' in para.keys():
    filemode = int(para['filemode'][0])
if "allpairinput" in para.keys():
    print("Use allpair.arc as input")
    filename = "allpair.arc"
    filemode = 2

if 'readallmin' in para.keys():
    Lallmin = int(para['readallmin'][0])

if 'LGibbs' in para.keys():
    LGibbs = int(para['LGibbs'][0])

if "multi" in para.keys():
    if 'surface' in para.keys():
        test = AllPair_surfacemode(LConly=LConly)
        test.multi = True
        test.LGibbs = LGibbs
        print('Pathanalyze start at %s' % (datetime.now()))
        print('Mode: surface-Multi   inputfilemode: %d  readmin: %d\n' %
              (filemode, Lallmin))
    else:
        # test = AllPair()
        print('Pathanalyze start at %s' % (datetime.now()))
        print('Mode: organic-Multi(NULL!)   inputfilemode: %d  readmin: %d\n' %
              (filemode, Lallmin))

else:
    if 'surface' in para.keys():
        test = AllPair_surfacemode(LConly=LConly)
        print('Pathanalyze start at %s' % (datetime.now()))
        print('Mode: surface   inputfilemode: %d  readmin: %d\n' %
              (filemode, Lallmin))
    else:
        # test = AllPair()
        print('Pathanalyze start at %s' % (datetime.now()))
        print('Mode: organic(NULL!)   inputfilemode: %d  readmin: %d\n' %
              (filemode, Lallmin))

if "db" in para.keys():
    # storage5 /home11/shiyf/program
    dbLoader = LaspDBPathSmpl("xxx", "xxxx", "xxx", "xxxxxxxx")

    rootDir = os.getcwd()
    globalVar.ReadLaspIn()
    print(globalVar.globalDict)

    if "dropDB" in para.keys():
        print("Going to drop current directory data in data base")
        dbLoader.PathSmplDrop(rootDir)
        sys.exit(0)
    if "multi" not in para.keys():
        AnalyzerCalOrDown(dbLoader, test, rootDir)
    else:
        # if multi in para:
        with os.popen("ls -d *.*/") as fp:
            dirs = [s.rstrip() for s in fp.readlines()]
        order = {"1": 1, "3": 2, "2": 3}
        dirs.sort(key=lambda x: order[re.match(r'^(\d+)', x).group(0)])
        alltest = []
        for direc in dirs:
            os.chdir(direc)
            if "surface" in para.keys():
                testTmp = AllPair_surfacemode(LConly=LConly)
            else:
                testTmp = None
            AnalyzerCalOrDown(dbLoader, testTmp, os.path.join(rootDir, direc))
            alltest.append(testTmp)
            os.chdir(rootDir)
        test.JoinAllPair(alltest)

else:
    if "multi" not in para.keys():
        test.readfile(filemode=filemode, Lallmininput=Lallmin, LGibbs=LGibbs)
        test.GibbsCorrectAll()
        test.GenPairdict_byname()
    else:
        rootDir = os.getcwd()
        with os.popen("ls -d *.*/") as fp:
            dirs = [s.rstrip() for s in fp.readlines()]
        order = {"1": 1, "3": 2, "2": 3}
        dirs.sort(key=lambda x: order[re.match(r'^(\d+)', x).group(0)])
        alltest = []
        for direc in dirs:
            os.chdir(direc)
            if "surface" in para.keys():
                testTmp = AllPair_surfacemode(LConly=LConly)
            else:
                testTmp = None
            testTmp.readfile(filemode=filemode,
                             Lallmininput=Lallmin,
                             LGibbs=LGibbs)
            testTmp.GibbsCorrectAll()
            testTmp.GenPairdict_byname()
            alltest.append(testTmp)
            os.chdir(rootDir)
        test.JoinAllPair(alltest)

# stringflag = True
stringflag = 0
lpts = 1
print("read file finish at  %s" % datetime.now())

if 'pts' in para.keys():
    lpts = int(para['pts'][0])

if 'stringflag' in para.keys():
    stringflag = int(para['stringflag'][0])

if 'all' in para.keys():
    print('job: out all minstr')
    test.OutInfo_script()

if 'ly' in para.keys():
    if len(para['ly']) == 1:
        sub = "sub.arc"
        substrate = allstr()
        substrate.readfile(sub)
        inputfile = para['ly'][0]
        name1, name2 = test.GetTarget(inputfile, substrate)
    if len(para['ly']) == 2:
        name1 = test.allnameid[int(para['ly'][0])]
        name2 = test.allnameid[int(para['ly'][1])]

    if "multi" in para.keys():
        name1 = test.multiEcfpNameDict[name1]
        name2 = test.multiEcfpNameDict[name2]

    print('test pair ECFP calculated. %s' % (datetime.now()))
    print('job: Yen try to link %s   %s \n' % (name1, name2))
    test.findpath_Yen(name1, name2, printstring=stringflag)
    for path in test.allpath[:10]:
        print(path.barrierlist)

if 'lmk' in para.keys():
    if len(para['lmk']) == 1:
        sub = "sub.arc"
        substrate = allstr()
        substrate.readfile(sub)
        inputfile = para['lmk'][0]
        name1, name2 = test.GetTarget(inputfile, substrate)
    if len(para['lmk']) == 2:
        name1 = test.allnameid[int(para['lmk'][0])]
        name2 = test.allnameid[int(para['lmk'][1])]

    if "multi" in para.keys():
        name1 = test.multiEcfpNameDict[name1]
        name2 = test.multiEcfpNameDict[name2]

    print('test pair ECFP calculated. %s' % (datetime.now()))
    print('job: Yen try to link %s   %s \n' % (name1, name2))
    test.findpath_MK(name1, name2, printstring=stringflag)


if 'mk' in para.keys():
    if len(para['mk']) == 1:
        inputfile = para['mk'][0]
        name1, name2 = test.GetTarget(inputfile)
    if len(para['mk']) == 2:
        name1 = test.allnameid[int(para['mk'][0])]
        name2 = test.allnameid[int(para['mk'][1])]

    if "multi" in para.keys():
        name1 = test.multiEcfpNameDict[name1]
        name2 = test.multiEcfpNameDict[name2]

    print('test pair ECFP calculated. %s' % (datetime.now()))
    test.findpath_MKMCXX(name1, name2, printstring=stringflag)

if 'mk-post' in para.keys():
    if len(para['mk-post']) == 1:
        inputfile = para['mk-post'][0]
        name1, name2 = test.GetTarget(inputfile)
    if len(para['mk-post']) == 2:
        name1 = test.allnameid[int(para['mk-post'][0])]
        name2 = test.allnameid[int(para['mk-post'][1])]

    if "multi" in para.keys():
        name1 = test.multiEcfpNameDict[name1]
        name2 = test.multiEcfpNameDict[name2]

    print('test pair ECFP calculated. %s' % (datetime.now()))
    test.findpath_MKMCXX_Post(name1, name2, printstring=stringflag)

if "kp" in para.keys():
    print('Select Key Pair %s' % (datetime.now()))
    test.PrintKeyPair(para["kp"][0])

if "kptarget" in para.keys():
    print('Select Key Pair %s' % (datetime.now()))
    test.PrintKeyPairTarget(para["kptarget"][0])

if 'l' in para.keys():
    if len(para['l']) == 1:
        inputfile = para['l'][0]
        name1, name2 = test.GetTarget(inputfile)
    if len(para['l']) == 2:
        name1 = test.allnameid[int(para['l'][0])]
        name2 = test.allnameid[int(para['l'][1])]
    print('test pair ECFP calculated. %s' % (datetime.now()))
    print('job: try to link %s   %s \n' % (name1, name2))
    test.findpath(name1, name2, printstring=stringflag)

if 'artifact' in para.keys():
    name1 = para['artifact'][0]
    name2 = para['artifact'][1]
    print('try to link %s %s' % (name1, name2))
    test.findpath(name1, name2, printstring=stringflag)

if 'targetpair' in para.keys():
    if len(para['targetpair']) == 1:
        inputfile = para['targetpair'][0]
        name1, name2 = test.GetTarget(inputfile)
    if len(para['targetpair']) == 2:
        name1 = test.allnameid[int(para['targetpair'][0])]
        name2 = test.allnameid[int(para['targetpair'][1])]

    test.PrintPair(name1, name2, lpts)

if 'tplist' in para.keys():
    test.PrintPair_list(para['tplist'][0], lpts)

if 'targetmin' in para.keys():
    if len(para['targetmin']) == 1:
        inputfile = para['targetmin'][0]
        test.OutTargetMin(inputfile)

if 'p' in para.keys():
    try:
        id = int(para['p'][0])
        test.OutPairdict_byname(test.allnameid[id])

    except Exception:
        Str = para['p'][0]
        name = test.GetSingleTarget(Str)
        test.OutPairdict_byname(name)

if 's' in para.keys():
    try:
        id = int(para['s'][0])
        test.OutAllPathfromSpestr(test.allnameid[id])

    except Exception:
        Str = para['s'][0]
        name = test.GetSingleTarget(Str)
        test.OutAllPathfromSpestr(name)

if 'allpair' in para.keys():
    print('job: out all minpair')
    test.Outallpair(pts=lpts)

if 'trainsetpair' in para.keys():

    import random
    size = int(para['trainsetpair'][0])
    print('job: make train set of pair. size = ' + str(size))
    _tmp = allstr()
    for pair in test.minpair.values():
        _tmp.append(pair.TSstr)
    print("allpair.arc contain:\n" + str(len(_tmp)) + " TS")
    if len(_tmp) > size:
        random.shuffle(_tmp)
        _tmp = _tmp[0:int(0.9 * size)]
        print("shffule " + str(len(_tmp)) + " TS")

    size = min(int((size - len(_tmp)) / 2), len(test.minpair.values()))
    for pair in random.sample(test.minpair.values(), size):
        _tmp.append(pair[0])
    for pair in random.sample(test.minpair.values(), size):
        _tmp.append(pair[1])
    print(str(size) + " IS\n" + str(size) + " FS")
    _tmp.printall('trainset.arc')
if "dotPlot" in para.keys():
    test.OutDotRsGraph()

print('Pathanalyze finish at %s' % (datetime.now()))
