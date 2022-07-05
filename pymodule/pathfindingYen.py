import mpmath
import os
from analyze_surfacemode import AllPair as AllPair_surfacemode
from analyze_surfacemode import GibbsCorrDict
from ECFP import AllStr as AllStr_ECFP
from datetime import datetime
import re
import numpy as np
import pprint
from multiprocessing import Pool, TimeoutError
from mkSolver import Spec, Reac, KbT, KbTh, mkBaseClass, mkSinglePoint

outputLevel = 1


def PathanalyzeYen(startName, sampledNormSminame, ballanceName=""):
    Lallmin = 2
    filemode = 1
    LGibbs = 500

    analyzer = AllPair_surfacemode(LConly=False)
    print('Pathanalyze start at %s' % (datetime.now()))
    print('Mode: surface   inputfilemode: %d  readmin: %d\n' %
          (filemode, Lallmin))
    analyzer.readfile(filemode=filemode, Lallmininput=Lallmin,
                      LGibbs=LGibbs)  # poolcount += 96*2
    analyzer.GibbsCorrectAll()
    analyzer.GenPairdict_byname()

    nameList = []
    for name2 in analyzer.allname:

        struc = analyzer.allstrmin[name2]
        fragments = [
            re.sub(r"-.*ads", "", frag) for frag in struc.sminame.split('.')
        ]
        # for HCOO
        # if "HCOO" not in fragments:
        #     continue

        _normSmiName2 = fragments + ballanceName.split(".")
        normSminame2 = ".".join(sorted(_normSmiName2))
        if normSminame2 in sampledNormSminame:
            print(normSminame2 + " sampled!")
            continue
        else:
            nameList.append(name2)
            print(normSminame2)
    with Pool(processes=24) as pool:
        jobs = []
        for name2 in nameList:
            j = pool.apply_async(analyzer.findpath_Yen,
                                 args=(startName, name2, False, False, 20))
            jobs.append([name2, j])
        try:
            allResult = [(n, j.get(timeout=60)) for n, j in jobs]
        except TimeoutError:
            allResult = [(n, j.get()) for n, j in jobs if j.ready()]

    # allResult = []
    # for name2 in nameList:
    #     print(name2, analyzer.allstrmin[name2].sminame)
    #     path = analyzer.findpath_Yen(startName, name2, False, False, 20)
    #     allResult.append((name2, path))

    # poolcount += 24
    allPath = []
    for n, r in allResult:
        if len(r) == 0:
            continue
        allPath.append([n, r[0]])
    if len(allPath) == 0:
        return None, None
    allPath.sort(key=lambda x: x[1].score)
    if outputLevel > 1:
        print(allPath)

    if outputLevel > 0:
        for name2, path in allPath[:10]:
            print("%.3e\t%s" % (path.score, analyzer.allstrmin[name2].sminame))

    nextSpe = allPath[0][0]
    print("nextSpe: %s" % analyzer.allstrmin[nextSpe].sminame)
    return nextSpe, analyzer.allstrmin[nextSpe]


if __name__ == "__main__":
    outputLevel = 1
    start = AllStr_ECFP()
    start.readfile("start.arc")
    start.GetAllSminame_FromEcfp(numproc=8, colorflag=0)
    start.GetAllECFPname()
    os.chdir("analyze")
    PathanalyzeYen(start[0].ECFPname, [start[0].ECFPname])
