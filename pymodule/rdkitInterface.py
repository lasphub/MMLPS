import numpy as np
# from rdkit import Chem
# from rdkit import RDLogger
# from bondOrder import Bond
# from multiprocessing import Pool
# from Font import outBlue, outSky, endMark
# import ctypes
# from ctypes import pointer
# from PeriodicTable import Elemass
from PeriodicTable import Eletable, Elemass
from copy import deepcopy

# RDLogger.DisableLog('rdApp.*')

readableNameDict = {
    "[CH2:1]=[O:2]": "CO",
    "[C:1](=[O:2])=[O:3]": "CO2",
    "[H:1][H:2]": "H2",
    "[HH:1]": "H",
    "[OH2:1]": "O",
    "[OH2:2]": "OH",
    "[OH2:3]": "H2O",
    "[CH:2](=[O:3])[OH:4]": "HCOO/COOH",
    "[CH:2]([OH:3])=[O:4]": "HCOO/COOH",
    "[CH:3](=[O:4])[OH:5]": "HCOOH",
    "[CH:3]([OH:4])=[O:5]": "HCOOH",
    "[CH2:3]([OH:4])[OH:5]": "H2CO/C(OH)2",
    "[CH2:4]([OH:5])[OH:6]": "H2COOH/HC(OH)2",
    "[CH2:4]([OH:5])[OH:6]-C1": "H2COOH/HC(OH)2",
    "[CH2:5]([OH:6])[OH:7]": "H2C(OH)2",
    "[CH2:2]=[O:3]": "CH=O",
    "[CH3:2][OH:3]": "COH",
    "[CH3:3][OH:4]": "HCOH",
    "[CH2:3]=[O:4]": "HCHO",
    "[CH3:4][OH:5]": "CH3O/CH2OH",
    "[CH3:5][OH:6]": "CH3OH",
    '[CH4:1]': "C",
    "[CH4:3]": "CH2",
    "[CH4:4]": "CH3",
    "[CH4:5]": "CH4",
}


def GenMol4Group(substr, bmx):
    allbond = []

    for i in range(len(substr)):
        for j in range(i + 1, len(substr)):
            if bmx[substr[i].id][substr[j].id] > 0:
                # if substr[i].elesymbol == "H" or substr[j].elesymbol == "H" : continue
                allbond.append([i, j, bmx[substr[i].id][substr[j].id]])
                # print(substr[i].elesymbol, substr[j].elesymbol,
                #       bmx[substr[i].id][substr[j].id])

    lines = []
    lines.append('%d  %d\n' % (0, 0))
    lines.append('     RDKit \n')  # %15.9f\n'%(str.energy))
    lines.append('\n')
    lines.append('%3d%3d  0  0  0  0  0  0  0  0999 V2000\n' %
                 (len(substr), len(allbond)))
    for i, atom in enumerate(substr):
        # if atom.elesymbol == "H":
        #     continue
        lines.append(
            '%10.4f%10.4f%10.4f%3s 0  0  0  0  0  0  0  0  0%3d  0  0\n' %
            (atom.xyz[0], atom.xyz[1], atom.xyz[2], atom.elesymbol, i + 1))
    for bond in allbond:
        lines.append('%3d%3d%3d  0\n' % (bond[0] + 1, bond[1] + 1, bond[2]))
    lines.append('M  END\n')
    lines.append('$$$$\n')
    # if len(lines) == 6: return None
    # if len(substr) > 1:
    #     # print("".join(lines))
    return "".join(lines)


# def sminame(substr, bmatrix, lat, bondneed=[], modeflag=1, surface=[]):
def RdkitSminame(substr, bmatrix, bondneed=[], modeflag=1, surface=[]):
    # judge whether a well molecular
    surfaceflag = '-'

    if modeflag == 1:
        pass
        # iza = [atom.ele for atom in substr]
        # sna = len(substr)**2
        # latinv = np.linalg.inv(lat)

        # bmx = []
        # for atom in substr:
        #     bmx.append([bmatrix[atom.id][jatm.id] for jatm in substr])
        # bmx = np.array(bmx).reshape(sna)

        # radical = [0 for atom in substr]
        # na = len(substr)
        # fa = reduce(lambda a, b: a + b,
        #             [list(np.matmul(atom.xyz, latinv)) for atom in substr])

        # c_na = pointer(ctypes.c_int(na))
        # c_iza = pointer((ctypes.c_int * na)(*iza))
        # c_rad = pointer((ctypes.c_int * na)(*radical))
        # c_bmx = pointer((ctypes.c_int * (sna))(*bmx))
        # c_fa = pointer((ctypes.c_double * (3 * na))(*fa))

        # c_lat = pointer(
        #     (ctypes.c_double * 9)(*reduce(lambda a, b: a + b, lat)))
        # # print 'test1'
        # program = ctypes.cdll.LoadLibrary(
        #     '/home10/kpl/pymodule/script/Lib_bondmatrix/bondmatrix.so')
        # # print 'test2'
        # program.hradicalorder_(c_na, c_iza, c_bmx, c_rad, c_fa, c_lat)
        # # print np.array(c_rad.contents)
        # if np.array(c_rad.contents).max() > 0:
        #     flag = False
        #     radflag = '-'
        #     for iatom, irad in enumerate(np.array(c_rad.contents)):
        #         if irad > 0:
        #             radflag = radflag + '%s%d' % (Eletable[iza[iatom] - 1],
        #                                           c_rad.contents[iatom])
        # else:
        #     flag = True

    if modeflag == 2:
        radflag = '-'
        flag = True
        for atom in substr:
            if bondneed[atom.id] != 0:
                flag = False
                radflag = radflag + '%s%d' % (Eletable[atom.ele - 1],
                                              bondneed[atom.id])

        for atom in substr:
            if surface[atom.id] == 1:
                surfaceflag = surfaceflag + '%sads' % (Eletable[atom.ele - 1])

    block = GenMol4Group(substr, bmatrix)
    mol = Chem.MolFromMolBlock(block)
    # m3 = Chem.RemoveHs(mol)
    if mol is not None:
        rdkitsmi = Chem.MolToSmiles(mol,
                                    allBondsExplicit=False,
                                    allHsExplicit=True)
        molecular = readableNameDict.get(rdkitsmi, rdkitsmi)

        if __name__ == "__main__" and rdkitsmi != "[Cu:1]" and rdkitsmi != "[Zn:1]":
            print((rdkitsmi, molecular))
    else:
        molecular = "UNKNOWN"
    # if molecular == '[H][H]':
    #     molecular = 'Hydrogen_Gas'
    # elif molecular == 'C(=O)=O':
    #     molecular = 'Carbon_Dioxide'
    # elif molecular == 'C(=O)O':
    #     if len(substr) == 5:
    #         molecular = 'Formic_Acid'
    #     elif len(substr) == 4:
    #         molecular = 'Formate/Carboxyl'
    # elif molecular == 'C':
    #     if len(substr) == 5:
    #         molecular = 'Methane'
    #         flag = True
    #     elif len(substr) == 4:
    #         molecular = 'Methyl'
    #         # flag = False
    # elif molecular == 'O':
    #     if len(substr) == 3:
    #         molecular = 'Water'
    #         flag = True
    #     elif len(substr) == 2:
    #         molecular = 'Hydroxy'
    #         # flag = False
    # elif molecular == 'C#[O]':
    #     molecular = 'Carbon_Monoxide'
    #     flag = True
    # elif molecular == 'C(=O)(O)O' and len(substr) == 6:
    #     molecular = 'Carbonic-Acid'
    # elif molecular == 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O':
    #     molecular = 'alpha_D_glucopyranose'
    # elif molecular == '[C@@H]1([C@@H]([C@H]([C@@H]([C@@H](CO)O1)O)O)O)O':
    #     molecular = 'beta_D_glucopyranose'
    # elif molecular == 'C(=O)[C@@H]([C@H]([C@@H]([C@@H](CO)O)O)O)O':
    #     molecular = 'P1 : D-glucose'
    # elif molecular == 'C(C(=O)[C@H]([C@@H]([C@@H](CO)O)O)O)O':
    #     molecular = 'D-fructose'
    # elif '[C@H]12[C@@H]([C@H]([C@@H]([C@@H](CO2)O1)O)O)O' in molecular:
    #     # elif molecular == '[C@H]12[C@@H]([C@H]([C@@H]([C@@H](CO2)O1)O)O)O + Water'        :
    #     molecular = 'P2 : Levoglucosan'
    # elif '[C@@H]1([C@@H]([C@@H]2[C@@H]([C@@H](CO2)O1)O)O)O' in molecular:
    #     molecular = 'P6 : Levoglucosan-2'
    # elif molecular == '[C@@H]1([C@@H]([C@H]([C@@H]([C@@H](CO)O)O1)O)O)O':
    #     molecular = 'P3 :  beta-D-glucofuranose '
    # elif molecular == '[C@@H]1([C@@H]([C@H]([C@@H]([C@@H](CO1)O)O)O)O)O':
    #     molecular = 'P4 : beta-D-glucoseptanose'
    # elif 'C(=O)[C@@H]1[C@H]([C@@H]([C@@H](CO)O1)O)O' in molecular:
    #     molecular = 'P5 : 2,5-anhydro-D-mnnose'

    # elif 'C(=O)[C@H]1[C@H]([C@@H]([C@@H](CO)O1)O)O' in molecular:
    #     molecular = 'P5 : 2,5-anhydro-D-mnnose-2'
    # # elif molecular == '' :
    # #    molecular =
    # elif molecular == 'CC':
    #     molecular = "Ethane"
    #     if len(substr) == 7:
    #         molecular = 'Ethyl'
    # elif molecular == 'O=O':
    #     molecular = 'Oxygen_Gas'
    # elif molecular == 'N#N':
    #     molecular = 'Nitrogen_Gas'
    # elif molecular == 'C=C':
    #     molecular = "Ethene"
    #     if len(substr) == 5:
    #         molecular = 'Ethenyl'
    # elif molecular == 'CO':
    #     if len(substr) == 6:
    #         molecular = "Methanol"
    # elif molecular == 'CCO' or molecular == 'C(C)O':
    #     molecular = "Ethanol"
    # elif molecular == 'CCC' or molecular == 'C(C)C':
    #     molecular = "Propane"
    #     if len(substr) == 10:
    #         molecular = 'Propyl'
    # elif molecular == 'C=O':
    #     if len(substr) == 2:
    #         molecular = 'Carbon_Monoxide'
    #         flag = True
    #     elif len(substr) == 4:
    #         molecular = 'Formaldehyde'
    #         flag = True
    # elif molecular == '[H]':
    #     flag = True
    # elif molecular == '[OH3]':
    #     flag = True

    if not flag:
        molecular = molecular + radflag
    if surfaceflag != '-':
        molecular = molecular + surfaceflag

    for i, atom in enumerate(substr):
        if atom.ele > 18:
            printmolecular = ''
            return printmolecular, flag
    # print(molecular)
    return molecular, flag


def calAllName(data):
    group, bondmatrix, bondneed, flag, surface = data
    return [
        RdkitSminame(sub, bondmatrix, bondneed, flag, surface) for sub in group
    ]


# def calAllpureName(data):
#     group, bondmatrix, lat, bondneed, flag, surface = data
#     return [
#         puresminame(sub, bondmatrix, lat, bondneed, flag, surface)
#         for sub in group
#     ]
#     # return [puresminame(group, bondmatrix, lat,bondneed,flag,surface)]


def calmass(str):
    allmass = np.array([Elemass[atm.ele - 1] for atm in str])
    return allmass.sum()


# def judgeReaction (strin1, strin2):


def singleFindName(strin):
    # str = deepcopy(strin)
    str = strin
    bmx = str.Bondmatrix()
    bmx2D = np.array(bmx).reshape(str.natm, str.natm)
    str.bmx1D = bmx
    group = str.Segmolecular(recal=False)
    substr = [[] for i in np.unique(group)]
    for id, atom in enumerate(str.atom):
        atom.id = id
        substr[group[id] - 1].append(atom)
    substr = sorted(substr, key=lambda x: calmass(x), reverse=True)
    return [RdkitSminame(sub, bmx2D, str.lat) for sub in substr], bmx2D


def glueSegStr(allmol):
    outstr, cellflag = "", True
    p = 0
    for i, (name, flag) in enumerate(allmol):
        if flag:
            font = outSky
        else:
            font = outBlue
        if name == '':
            continue
        if p > 0:
            outstr += "+%s%s%s" % (font, name, endMark)
        elif p == 0:
            outstr += "%s%s%s" % (font, name, endMark)
            p = p + 1
        if not flag:
            cellflag = False
    return outstr.strip(), cellflag


def glueSegStr_pure(allmol):
    outstr, cellflag = "", True
    p = 0
    for i, (name, flag) in enumerate(allmol):
        if name == '':
            continue
        if p > 0:
            outstr += ".%s" % (name)
        elif p == 0:
            outstr += "%s" % (name)
            p = p + 1
        if not flag:
            cellflag = False
    return outstr.strip(), cellflag


if __name__ == "__main__":
    import allstr_new
    _tmp = allstr_new.allstr()
    # _tmp.readfile("/home10/shiyf/prog/laspy/Basic/test/test_MultStruc.arc")
    _tmp.readfile(
        "/home11/shiyf/5.ZnOH/0614_asop3/1.ZnOH111/../workspace-primslab/cu111/num-6/super_3xroot3/Zn-3/O-4/H-5/GM.arc"
    )
    _tmp.calAllFakeBondMaxtrix(numproc=6)
    _tmp.calAllSegmolecular(numproc=6)

    allgroup = []
    for str in _tmp:
        substr = [[] for i in np.unique(str.group)]

        for id, atom in enumerate(str.atom):
            atom.id = id
            substr[str.group[id] - 1].append(atom)
        substr = sorted(substr, key=lambda x: calmass(x), reverse=True)
        allgroup.append((substr, str.bmx2D, str.bondneed, 2, str.surfaceatom))

    for i in range(0, len(allgroup)):
        print(glueSegStr_pure(calAllName(allgroup[i])))
