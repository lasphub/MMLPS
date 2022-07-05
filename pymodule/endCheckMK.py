import mpmath

from analyze_surfacemode import AllPair as AllPair_surfacemode
from analyze_surfacemode import GibbsCorrDict
from ECFP import AllStr as AllStr_ECFP
from datetime import datetime
import re
import numpy as np
import pprint
from multiprocessing import Pool, TimeoutError
from mkSolver import Spec, Reac, KbT, KbTh, mkBaseClass, mkSinglePoint
import os

allowedGas = ["co2", "co", "ch3oh", "h2", "h2o"]
outputLevel = 2


def Smi2SpeName(sminame):
    n = sminame.lower()
    if "ads" in n:
        return n.split("-")[0] + "*"
    elif n not in allowedGas:
        if outputLevel > 2:
            print("WARNING: %s is treated as adsorbates" % n)
        return n + "*"
    else:
        return n


def GenMKModel(analyzer, mkPair):

    # reacNameExist = []
    species = []
    reaction = []
    for pair in mkPair:
        bf = pair.TSstr.energy - analyzer.allstrmin[pair[0].name].energy
        bb = pair.TSstr.energy - analyzer.allstrmin[pair[1].name].energy

        reactant = analyzer.allstrmin[pair[0].name].sminame.split(".")
        product = analyzer.allstrmin[pair[1].name].sminame.split(".")
        tstate = pair.TSstr.sminame.split(".")
        if "HCOOH-Oads" in reactant and "HCOOH" in tstate:
            bf -= GibbsCorrDict["HCOOH"]
        elif "HCOOH" in reactant and "HCOOH-Oads" in tstate:
            bf += GibbsCorrDict["HCOOH"]
        elif "HCOOH-Oads" in product and "HCOOH" in tstate:
            bb -= GibbsCorrDict["HCOOH"]
        elif "HCOOH" in product and "HCOOH-Oads" in tstate:
            bb += GibbsCorrDict["HCOOH"]

        if "HCHO-CadsOads" in reactant and "HCHO" in tstate:
            bf -= GibbsCorrDict["HCHO"]
        elif "HCHO" in reactant and "HCHO-CadsOads" in tstate:
            bf += GibbsCorrDict["HCHO"]
        elif "HCHO-CadsOads" in product and "HCHO" in tstate:
            bb -= GibbsCorrDict["HCHO"]
        elif "HCHO" in product and "HCHO-CadsOads" in tstate:
            bb += GibbsCorrDict["HCHO"]

        bystander = []
        # pick out bystanders
        pair.bystander = bystander
        while True:
            flag = False
            for r in reactant:
                if r in product:
                    flag = True
                    break
            if flag is False:
                break
            reactant.remove(r)
            product.remove(r)
            bystander.append(r)
        # print(".".join(reactant), ".".join(product), ".".join(bystander))
        if len(reactant) == 0 or len(product) == 0:
            pair.bf = 0
            pair.bb = 0
            pair.reactant = []
            pair.product = []
            pair.bystander = []
            continue

        reactant = [Smi2SpeName(s) for s in reactant]
        product = [Smi2SpeName(s) for s in product]
        reactSite = "".join(reactant).count("*")
        prodSite = "".join(product).count("*")
        if reactSite > prodSite:
            product.extend(["*"] * (reactSite - prodSite))
        elif prodSite > reactSite:
            reactant.extend(["*"] * (prodSite - reactSite))

        pair.bf = bf
        pair.bb = bb
        pair.reactant = reactant
        pair.product = product
        pair.bystander = bystander
        pair.name = ".".join(sorted(reactant + product))

    mkPair_sim = {}
    for pair in mkPair:
        if pair.name not in mkPair_sim:
            mkPair_sim[pair.name] = pair
        else:
            if outputLevel > 2:
                print("WARNING: dummy reaction: %s: %s--%s && %s--%s" % (
                    pair.name,
                    analyzer.allstrmin[pair[0].name].sminame,
                    analyzer.allstrmin[pair[1].name].sminame,
                    analyzer.allstrmin[mkPair_sim[pair.name][0].name].sminame,
                    analyzer.allstrmin[mkPair_sim[pair.name][1].name].sminame,
                ))

    mkBaseClass.energyType = "barrier"
    species = set()
    for rn in mkPair_sim.keys():
        species.update(set(rn.split(".")))
    # speName = {s: s.lower().split("-")[0] for s in species}

    species = {n: Spec(n) for n in species}
    add_species = {
        "*": Spec("*"),
        "co2": Spec("co2"),
        "co": Spec("co"),
        "co*": Spec("co*"),
        "ch3oh": Spec("ch3oh"),
        "ch3oh*": Spec("ch3oh*"),
        "h2": Spec("h2"),
        "h2o": Spec("h2o"),
        "h2o*": Spec("h2o*"),
    }
    species.update(add_species)

    adsAddBarrier = 0.3
    reaction = [
        # Reac([species["co2"], species["*"]], [species["co2*"]], bf=0, bb=0),
        # Reac([species["co"], species["*"]], [species["co*"]], bf=0.19, bb=0),
        Reac([species["co"], species["*"]], [species["co*"]],
             bf=0.25 + adsAddBarrier,
             bb=0 + adsAddBarrier),
        Reac([species["h2o"], species["*"]], [species["h2o*"]],
             bf=0.40 + adsAddBarrier,
             bb=0 + adsAddBarrier),
        Reac([species["ch3oh"], species["*"]], [species["ch3oh*"]],
             bf=0.45 + adsAddBarrier,
             bb=0 + adsAddBarrier),
    ]

    for p in mkPair_sim.values():
        if len(p.reactant) == 0:
            continue
        reactant = [species[s] for s in p.reactant]
        product = [species[s] for s in p.product]
        reac = Reac(reactant, product, bf=p.bf, bb=p.bb)
        reac.sourcePair = p
        reaction.append(reac)

    for p in mkPair:
        for r in reaction:
            if hasattr(r, "sourcePair") and pair.name == r.sourcePair.name:
                p.reac = r
    return species, reaction


def RunMkModel(species, reaction):
    mkBaseClass.output = -1
    # gasPressure = {"h2": 40, "co": 10, "co2": 10, "ch3oh": 1, "h2o": 1}
    gasPressure = {"h2": 1, "co": 1, "co2": 1, "ch3oh": 1, "h2o": 1}
    temperature = 500

    # for k, v in species.items():
    #     if k in saveInitCGuess:
    #         v.c = saveInitCGuess[k]
    try:
        maxMkCycle = 3
        mkCycle = 0
        while mkCycle < maxMkCycle:
            task = mkSinglePoint.SingleMK(species, reaction, temperature,
                                          gasPressure)
            task.Run(odeTime=1e7, opt=True)
            if outputLevel > 2:
                pprint.pprint(task.result)

            if task.success:
                break
            for k, v in task.result[3].items():
                if v > 0:
                    species[k].c = v
                else:
                    species[k].c = -v / 1000
            mkCycle += 1
        for k, v in task.result[3].items():
            species[k].c = v
        task.success = True
        # saveInitCGuess[k] = v
    except TypeError:
        print("TypeError")
        task.success = False
        return task
    pprint.pprint(task.result)
    return task


def EndCheckMK(direcs, targetProd, jointAnaDirec):
    os.chdir(jointAnaDirec)
    direcs = [d.split("/")[-1] for d in direcs]
    for direc in direcs:
        fName = direc + "/allgoodsect.arc"
        if not os.path.isfile(fName):
            print("%s not exist" % fName)
            return False
    Lallmin = 2
    filemode = 1
    LGibbs = 500

    # analyzer = AllPair_surfacemode(LConly=False)
    # print('Pathanalyze start at %s' % (datetime.now()))
    # print('Mode: surface   inputfilemode: %d  readmin: %d\n' %
    #       (filemode, Lallmin))
    # analyzer.readfile(filemode=filemode, Lallmininput=Lallmin,
    #                   LGibbs=LGibbs)  # poolcount += 96*2
    # analyzer.GibbsCorrectAll()
    # analyzer.GenPairdict_byname()

    analyzer = AllPair_surfacemode(LConly=False)
    analyzer.multi = True
    analyzer.LGibbs = LGibbs
    print('Pathanalyze start at %s' % (datetime.now()))
    print('Mode: surface-Multi   inputfilemode: %d  readmin: %d\n' %
          (filemode, Lallmin))

    alltest = []
    for direc in direcs:
        os.chdir(direc)
        testTmp = AllPair_surfacemode(LConly=False)
        testTmp.readfile(filemode=filemode,
                         Lallmininput=Lallmin,
                         LGibbs=LGibbs)
        testTmp.GibbsCorrectAll()
        testTmp.GenPairdict_byname()
        alltest.append(testTmp)
        os.chdir(jointAnaDirec)
    analyzer.JoinAllPair(alltest)

    nameDict = {}
    for name2 in analyzer.allstrmin.keys():

        struc = analyzer.allstrmin[name2]
        normSminame2 = ".".join(
            sorted([
                re.sub(r"-.*ads", "", frag)
                for frag in struc.sminame.split('.')
            ]))
        nameDict[normSminame2] = name2
    # print(nameDict)
    targetProdName = []
    for reacProdPair in targetProd:
        normSminame1 = ".".join(
            sorted([
                re.sub(r"-.*ads", "", frag)
                for frag in reacProdPair[0].split('.')
            ]))
        normSminame2 = ".".join(
            sorted([
                re.sub(r"-.*ads", "", frag)
                for frag in reacProdPair[1].split('.')
            ]))
        # print(normSminame1, normSminame2)
        if normSminame1 not in nameDict:
            print("%s not in nameDict" % normSminame1)
            return False
        elif normSminame2 not in nameDict:
            print("%s not in nameDict" % normSminame2)
            return False
        targetProdName.append((nameDict[normSminame1], nameDict[normSminame2]))

    for reacProdPair in targetProdName:
        allPath = analyzer.findpath_Yen(reacProdPair[0], reacProdPair[1],
                                        False, False, 10)
        if len(allPath) != 0:
            # analyzer = _analyzer
            mkPair = list()
            for path in allPath:
                for pair in path:
                    if pair not in mkPair:
                        mkPair.append(pair)
            species, reaction = GenMKModel(analyzer, mkPair)
            task = RunMkModel(species, reaction)
            # for path in analyzer.allpath:
            #     print(path)
            print(analyzer.allpath[0])
            print("CH3OH rate: %e" % task.result[2]["ch3oh"])

            if task.success:
                return float(task.result[2]["ch3oh"])
            else:
                print("mk Fail")
                return 0
        else:
            print("No Path!")
            return False


if __name__ == "__main__":
    outputLevel = 3
    # direcs, targetReacProdList, jointAnaDirec
    EndCheckMK(["1.CO2-HCOOH", "3.HCHO-CH3OH", "2.HCOOH-HCHO"],
               [("CO2.H2.H2.H2", "CH3OH.H2O")],
               "/home10/shiyf/3.CO2Hydro/1115/3.211-0Zn")
    # EndCheckMK(["1.CO2+H2", "3.HCHO+H2", "2.HCOOH+H2"],
    #            [("CO2.H2.H2.H2", "CH3OH.H2O")],
    #            "/home10/shiyf/3.CO2Hydro/0320_test/joint")