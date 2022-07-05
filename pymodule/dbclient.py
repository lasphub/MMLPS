#! /home10/shiyf/bin/miniconda3/bin/python
import pymongo
# from bson.objectid import ObjectId
import os
from pprint import pprint
import argparse
from allstr_new import allstr
from analyze_surfacemode import Pair
from structure_new import Str as Structure
from ECFP import hashlist
import hashlib
import PeriodicTable

# def print_json(data):
#     print(json.dumps(data, sort_keys=True, indent=4, separators=(', ', ': '), ensure_ascii=False))


def GenArcString(strucList, head=False):
    """Gen Arc file without"""
    lines = []
    # if head is True:
    #     lines.append('!BIOSYM archive 2\n')
    #     lines.append('PBC=ON\n')
    for istr, struc in enumerate(strucList):
        istr = istr + 1
        if hasattr(struc, 'label'):
            label = struc.label
        else:
            label = ''
        lines.append('     Energy     %8d    %8d     %12.6f\n' %
                     (istr, istr, struc.trueEnergy))
        lines.append('!DATE  %s\n' % label)
        lines.append('PBC %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n' %
                     (struc.abc[0], struc.abc[1], struc.abc[2], struc.abc[3],
                      struc.abc[4], struc.abc[5]))
        # f.write('PBC %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n'%(100,100,100,90,90,90))
        for i, atom in enumerate(struc.atom):
            lines.append(
                '%-3s %15.9f %15.9f %15.9f CORE %4d %2s %2s %8.4f %4d\n' %
                (atom.elesymbol, atom.xyz[0], atom.xyz[1], atom.xyz[2], i + 1,
                 atom.elesymbol, atom.elesymbol, atom.charge, i + 1))
        lines.append('end\nend\n')
    return "".join(lines)


def GenStructureHash(struc, compose=None):
    if compose is None:
        compose = dict(
            zip([PeriodicTable.Eletable[k - 1] for k in struc.eleList],
                list(struc.natompe)))
    feature = []
    feature.append("".join("%d%s" % (v, k) for k, v in compose.items()))
    feature.extend(struc.abc)
    elelist = list(compose.keys())
    while len(elelist) != 0:
        for at in struc.atom:
            if at.elesymbol in elelist:
                feature.append(at.xyz[0])
                feature.append(at.xyz[1])
                feature.append(at.xyz[2])
                elelist.remove(at.elesymbol)
    strucHash = hashlist(feature, dim=10e10)
    return strucHash


# def GenPairHash(_pair):
#     compose = dict(
#         zip([PeriodicTable.Eletable[k - 1] for k in _pair[0].eleList],
#             list(_pair[0].natompe)))
#     feature = []
#     feature.append("".join("%d%s" % (v, k) for k, v in compose.items()))
#     feature.extend(_pair[0].abc)
#     for struc in [_pair[0], _pair.TSstr, _pair[1]]:
#         elelist = list(compose.keys())
#         while len(elelist) != 0:
#             for at in struc.atom:
#                 if at.elesymbol in elelist:
#                     feature.append(at.xyz[0])
#                     feature.append(at.xyz[1])
#                     feature.append(at.xyz[2])
#                     elelist.remove(at.elesymbol)
#     pairHash = hashlist(feature, dim=10e10)
#     return pairHash


class LaspDBBasic(object):
    def __init__(self, dbname, user, pwd, host, port=27017):
        dbUri = "mongodb://%s:%s@%s:%s/%s" % (user, pwd, host, port, dbname)
        print(dbUri)
        client = pymongo.MongoClient(dbUri)
        self.client = client

        self.infoCol = None
        self.fileCol = None
        self.infoID = None
        self.fileID = None

    def _DBUploadSingle(self, info, file):
        self.infoID = self.infoCol.insert_one(info).inserted_id
        self.fileID = self.fileCol.insert_one(file).inserted_id
        fUpRtn = self.fileCol.update_one({"_id": self.fileID},
                                         {"$set": {
                                             "info": self.infoID,
                                         }})
        iUpRtn = self.infoCol.update_one({"_id": self.infoID},
                                         {"$set": {
                                             "file": self.fileID,
                                         }})
        if fUpRtn.matched_count != 1 or iUpRtn.matched_count != 1:
            raise pymongo.errors.InvalidOperation


class LaspDBGlobalMin(LaspDBBasic):
    def __init__(self, dbname, user, pwd, host, port=27017):
        super(LaspDBGlobalMin, self).__init__(dbname,
                                              user,
                                              pwd,
                                              host,
                                              port=27017)

    def DBUploadStructure(self, struc, uploadType=0, **kwd):
        """
        Upload Global Minimum Structure.
        Necessary argument: category, name, author, exploreType
        """

        allowedKeys = [
            "category",
            "name",
            "description",
            "pot",
            "exploreType",
        ]
        if kwd["keys"]:
            for para in kwd["keys"]:
                p = para.split("=")
                if p[0] not in allowedKeys:
                    print("Err: " + p[0] + "  in not an allowed Key")
                else:
                    kwd[p[0]] = p[1]

        necessaryKwd = ["category", "name", "author", "exploreType"]
        flag = [False for key in necessaryKwd if key not in kwd]
        if flag is False:
            print(
                """Err: Missing Keywords. Provide them as "keys='name=methanol category=molecule'" or "name=methanol"
                    Following keyword is necessary: """ +
                "/".join * necessaryKwd)
            return

        db = self.client["GlobalMinDB"]
        self.fileCol = db["arcFile"]
        self.infoCol = db["strucInfo"]

        try:
            compose = dict(
                zip([PeriodicTable.Eletable[k - 1] for k in struc.eleList],
                    [float(f) for f in struc.natompe]))
            strucHash = GenStructureHash(struc, compose)
            info = {
                "category": kwd["category"],
                "name": kwd["name"],
                "description": kwd.get("description", ""),
                "author": kwd["author"],
                "pot": kwd.get("pot", ""),
                "exploreType": kwd["exploreType"],
                "compose": compose,
                "nnEnergy": struc.trueEnergy,
                "strucHash": strucHash,
                "dir": os.environ["HOSTNAME"] + ":" + os.getcwd()
            }
            file = {
                "file": GenArcString([struc]),
                "strucHash": strucHash,
            }
            self._DBUploadSingle(info, file)
            print("Success Upload Struc" +
                  "".join("%d%s" % (v, k)
                          for k, v in compose.items()) + "(%s)" % strucHash)
        except pymongo.errors.InvalidOperation as identifier:
            print("Err: InvalidOperation on DBUploadStructure")
            print(identifier)
            if self.fileID is not None:
                self.fileCol.delete_one({"_id": self.fileID})
            if self.infoID is not None:
                self.infoCol.delete_one({"_id": self.infoID})

    def DBDownloadStructure(self,
                            client,
                            dbname,
                            colname,
                            condition={},
                            file=False,
                            laspView=False,
                            **kwd):
        db = client[dbname]
        infoCol = db[colname]

        if kwd["keys"]:
            for para in kwd["keys"]:
                p = para.split("=")
                condition[p[0]] = p[1]
        print(condition)
        posts = infoCol.find(condition).limit(5)
        for post in posts:
            pprint(post, sort_dicts=False)
            if file:
                fname = "".join(
                    "%d%s" % (v, k)
                    for k, v in post["compose"].items()) + "-" + str(
                        post["strucHash"]) + ".arc"
                line = db["arcFile"].find_one({"_id": post["file"]})["file"]
                print("Write " + fname + " file")
                with open(fname, "w") as fp:
                    fp.write(line)


class LaspDBPathSmpl(LaspDBBasic):
    def __init__(self, dbname, user, pwd, host, port=27017):
        super(LaspDBPathSmpl, self).__init__(dbname,
                                             user,
                                             pwd,
                                             host,
                                             port=27017)
        self.allPair = None
        self.sumy = None

        self.namePairField = {
            "readNameIS",
            "readNameTS",
            "readNameFS",
            "ecfpNameMoleIS",
            "ecfpNameMoleTS",
            "ecfpNameMoleFS",
            "ecfpNameSurfIS",
            "ecfpNameSurfTS",
            "ecfpNameSurfFS",
            "ecfpNameSiteIS",
            "ecfpNameSiteTS",
            "ecfpNameSiteFS",
        }
        self.optionalPairField = [
            "reactEcfpNameMole",
            "reactEcfpNameSurf",
            "reactEcfpNameSite",
            "fpmatchLink",
            "fpmatchAtomMap",
            "reactCenterList0",
            "reactCenterList1",
            "siteMatchPair",
            "zpeIS",
            "zpeTS",
            "zpeFS",
            # "gasMoleIS",
            # "gasMoleTS",
            # "gasMoleFS",
        ]
        self.necessaryPairField = [
            "heatIF",
            "heatFI",
            "barrierIF",
            "barrierFI",
            "energyIS",
            "energyTS",
            "energyFS",
            "pathType",
            "strucHashIS",
            "strucHashTS",
            "strucHashFS"
            "fileIS",
            "fileTS",
            "fileFS",
            "sumyID",
        ]
        self.nameMinField = [
            "readNameLM",
            "ecfpNameMoleLM",
            "ecfpNameSurfLM",
            "ecfpNameSiteLM",
        ]
        self.optionalMinField = [
            "zpeLM",
            "gasMoleLM",
        ]
        self.necessaryMinField = [
            "energyLM",
            "pathType",
            "strucHash",
            "sumyID",
        ]

    # def __new__(self):

    def _GenpathSmplHash(self, dir):
        "fast read allgoodsect.arc/start.arc/, get an hash value and check if it is existed in Database"
        m = hashlib.md5()
        with open(os.path.join(dir, "allgoodsect.arc"), "rb") as f:
            while True:
                buf = f.read(2**20)
                if not buf:
                    break
                m.update(buf)
        with open(os.path.join(dir, "addmin.arc"), "rb") as f:
            while True:
                buf = f.read(2**20)
                if not buf:
                    break
                m.update(buf)
        return m.hexdigest()

    def PathSmplUploadCheckExist(self, dir):
        pathSmplHash = self._GenpathSmplHash(dir)
        self.pathSmplHash = pathSmplHash
        print(pathSmplHash)
        db = self.client["PathSmplDB"]
        sumyCol = db["Summary"]
        sumy = sumyCol.find_one({"pathSmplHash": pathSmplHash})
        print(sumy)
        if sumy is None:
            return False
        else:
            self.sumyCol = db["Summary"]
            self.sumy = sumy
            self.infoCol = db[sumy["ReactInfo"]]
            self.fileCol = db[sumy["ArcFile"]]
            return True

    def _PathSmplUploadPair(self, pair, sumyID):

        infoName = {}
        if hasattr(pair[0], "sminame") and hasattr(pair[1], "sminame"):
            infoName["readNameIS"] = pair[0].sminame
            infoName["readNameTS"] = pair.TSstr.sminame
            infoName["readNameFS"] = pair[1].sminame
        if hasattr(pair[0], "ECFPname") and hasattr(pair[0], "ECFPname"):
            infoName["ecfpNameMoleIS"] = pair[0].ECFPname
            infoName["ecfpNameMoleTS"] = pair.TSstr.ECFPname
            infoName["ecfpNameMoleFS"] = pair[1].ECFPname
        if hasattr(pair[0], "gasMole") and hasattr(pair[1], "gasMole"):
            infoName["gasMoleIS"] = pair[0].gasMole
            infoName["gasMoleTS"] = pair.TSstr.gasMole
            infoName["gasMoleFS"] = pair[1].gasMole
        if hasattr(pair.TSstr, "zpe"):
            infoName["zpeTS"] = pair.TSstr.zpe
        if hasattr(pair[0], "zpe"):
            infoName["zpeIS"] = pair[0].zpe
        if hasattr(pair[1], "zpe"):
            infoName["zpeFS"] = pair[1].zpe
        infoOpt = {}
        for attr in self.optionalPairField:
            if hasattr(pair, attr):
                infoOpt[attr] = pair.attr

        strucHashIS = GenStructureHash(pair[0])
        strucHashTS = GenStructureHash(pair.TSstr)
        strucHashFS = GenStructureHash(pair[1])
        info = {
            "energyIS": pair[0].trueEnergy,
            "energyTS": pair.TSstr.trueEnergy,
            "energyFS": pair[1].trueEnergy,
            "heatIF": pair[1].trueEnergy - pair[0].trueEnergy,
            "heatFI": pair[0].trueEnergy - pair[1].trueEnergy,
            "barrierIF": pair.TSstr.trueEnergy - pair[0].trueEnergy,
            "barrierFI": pair.TSstr.trueEnergy - pair[1].trueEnergy,
            "pathType": "pair",
            "strucHashIS": strucHashIS,
            "strucHashTS": strucHashTS,
            "strucHashFS": strucHashFS,
            "sumyID": sumyID,
        }
        fileIS = {
            "file": GenArcString([pair[0]]),
            "strucHash": strucHashIS,
        }
        fileTS = {
            "file": GenArcString([pair.TSstr]),
            "strucHash": strucHashTS,
        }
        fileFS = {
            "file": GenArcString([pair[1]]),
            "strucHash": strucHashFS,
        }
        info = {**infoName, **info, **infoOpt}
        self.infoID = self.infoCol.insert_one(info).inserted_id
        self.fileIDIS = self.fileCol.insert_one(fileIS).inserted_id
        self.fileIDTS = self.fileCol.insert_one(fileTS).inserted_id
        self.fileIDFS = self.fileCol.insert_one(fileFS).inserted_id
        fUpRtn = self.fileCol.update_one({"_id": self.fileIDIS},
                                         {"$set": {
                                             "info": self.infoID,
                                         }})
        fUpRtn = self.fileCol.update_one({"_id": self.fileIDTS},
                                         {"$set": {
                                             "info": self.infoID,
                                         }})
        fUpRtn = self.fileCol.update_one({"_id": self.fileIDFS},
                                         {"$set": {
                                             "info": self.infoID,
                                         }})
        iUpRtn = self.infoCol.update_one({"_id": self.infoID}, {
            "$set": {
                "fileIS": self.fileIDIS,
                "fileTS": self.fileIDTS,
                "fileFS": self.fileIDFS,
            }
        })
        if fUpRtn.matched_count != 1 or iUpRtn.matched_count != 1:
            raise pymongo.errors.InvalidOperation

        # print("Success Upload Pair" + infoName["readNameIS"] + " to " +
        #       infoName["readNameFS"] + " (%s)" % pair.name)

    def _PathSmplUploadMin(self, struc, sumyID):
        strucHash = GenStructureHash(struc)

        infoName = {}
        if hasattr(struc, "sminame"):
            infoName["readNameLM"] = struc.sminame
        if hasattr(struc, "ECFPname"):
            infoName["ecfpNameMoleLM"] = struc.ECFPname
        if hasattr(struc, "gasMole"):
            infoName["gasMoleLM"] = struc.gasMole
        if hasattr(struc, "zpe"):
            infoName["zpeLM"] = struc.zpe
        infoOpt = {}
        for attr in self.optionalMinField:
            if hasattr(struc, attr):
                infoOpt[attr] = struc.attr

        info = {
            "energyLM": struc.trueEnergy,
            "pathType": "min",
            "strucHashLM": strucHash,
            "sumyID": sumyID,
        }

        file = {
            "file": GenArcString([struc]),
            "strucHash": strucHash,
        }

        info = {**infoName, **info, **infoOpt}
        self.infoID = self.infoCol.insert_one(info).inserted_id
        self.fileID = self.fileCol.insert_one(file).inserted_id
        fUpRtn = self.fileCol.update_one({"_id": self.fileID},
                                         {"$set": {
                                             "info": self.infoID,
                                         }})
        iUpRtn = self.infoCol.update_one({"_id": self.infoID},
                                         {"$set": {
                                             "fileLM": self.fileID,
                                         }})
        if fUpRtn.matched_count != 1 or iUpRtn.matched_count != 1:
            raise pymongo.errors.InvalidOperation
        # print("Success Upload LocalMin " + infoName["readNameLM"] +
        #       " (%s)" % struc.name)

    def PathSmplUpload(self, analyzeAllPair, startSurf, startMole):
        self.allPair = analyzeAllPair

        db = self.client["PathSmplDB"]
        infoColName = "-".join(
            [startSurf, startMole, "ReactInfo", self.pathSmplHash])
        fileColName = "-".join(
            [startSurf, startMole, "ArcFile", self.pathSmplHash])
        sumyCol = db["Summary"]
        self.infoCol = db[infoColName]
        self.fileCol = db[fileColName]

        try:
            sumy = {
                "ReactInfo": infoColName,
                "ArcFile": fileColName,
                "startSurfReadName": startSurf,
                "startMoleReadName": startMole,
                "pathSmplHash": self.pathSmplHash,
                "author": os.environ["USER"],
                "dir": os.environ["HOSTNAME"] + ":" + os.getcwd(),
                "sourcedir": "?",
                "potentialType": "",
                "potentialPot": ""
                ""
            }
            sumyID = sumyCol.insert_one(sumy).inserted_id

            for _pair in self.allPair:
                self._PathSmplUploadPair(_pair, sumyID)

            for _struc in self.allPair.addmin:
                self._PathSmplUploadMin(_struc, sumyID)

        except pymongo.errors.InvalidOperation as identifier:
            print("Err: InvalidOperation on DBUploadStructure")
            print(identifier)
            if sumyID is not None:
                sumyCol.delete_one({"_id": sumyID})
            if self.fileCol is not None:
                self.fileCol.drop()
            if self.infoCol is not None:
                self.infoCol.drop()

    def _PathSmplDownloadStrucInit(self, info, suffix="", name="ecfpNameMole"):
        """Loading the info into the Structure Class Instance"""
        _struc = Structure()
        _struc.fileCol = self.fileCol
        _struc.name = info[name + suffix]
        attr = vars(_struc)
        for key in info:
            if key.rstrip(suffix) == "readName":
                attr["sminame"] = info[key]
            elif key.endswith(suffix):
                attr[key.rstrip(suffix)] = info[key]
        return _struc

    def PathSmplDownLoadInfo(self, analyzeAllPair, dir, name="ecfpNameMole"):
        """Dowload All Information from Database and stored them in analyzeAllPair"""
        if len(analyzeAllPair) != 0:
            print(
                "Err: analyzeAllPair is not empty and is not suitable for download info collection"
            )
            return
        if self.sumy is None:
            self.PathSmplUploadCheckExist(dir)
        self.allPair = analyzeAllPair
        allInfo = self.infoCol.find()

        for info in allInfo:
            if info["pathType"] == "pair":
                _strucIS = self._PathSmplDownloadStrucInit(info,
                                                           suffix="IS",
                                                           name=name)
                _strucTS = self._PathSmplDownloadStrucInit(info,
                                                           suffix="TS",
                                                           name=name)
                _strucFS = self._PathSmplDownloadStrucInit(info,
                                                           suffix="FS",
                                                           name=name)
                _pair = Pair(_strucIS, _strucFS, _strucTS.energy, _strucTS)
                pattr = vars(_pair)
                for key in info:
                    if not (key.endswith("IS") or key.endswith("TS")
                            or key.endswith("FS")):
                        pattr[key] = info[key]
                self.allPair.append(_pair)
            elif info["pathType"] == "min":
                _strucLM = self._PathSmplDownloadStrucInit(info,
                                                           suffix="LM",
                                                           name=name)
                self.allPair.addmin.append(_strucLM)
                print(_strucLM.sminame)

        for pair in self.allPair:
            for struc in pair:
                if struc.name not in self.allPair.allname.keys():
                    nameid = len(self.allPair.allname)
                    self.allPair.allname[struc.name] = nameid
                    self.allPair.allnameid[nameid] = struc.name
                    self.allPair.allstrbyname[struc.name] = [struc]
                else:
                    self.allPair.allstrbyname[struc.name].append(struc)

        for struc in self.allPair.addmin:
            if struc.name not in self.allPair.allname.keys():
                nameid = len(self.allPair.allname)
                self.allPair.allname[struc.name] = nameid
                self.allPair.allnameid[nameid] = struc.name
                self.allPair.allstrbyname[struc.name] = [struc]
            else:
                self.allPair.allstrbyname[struc.name].append(struc)

        self.allPair.Lallmin = 1

    def pathSmplDownLoadArc(self, struc):
        """Dowload Arcfile of giving ObjectID and return as string"""
        if self.fileCol is None and not hasattr(struc, "fileCol"):
            print("Err: infoCol or fileCol has not been set!")
            return
        if hasattr(struc, "fileCol"):
            return struc.fileCol.find_one({"_id": struc.file})["file"]
        else:
            return self.fileCol.find_one({"_id": struc.file})["file"]

    def PathSmplDrop(self, dir):
        if self.PathSmplUploadCheckExist(dir) is False:
            print("Data Not Exist. Exit.")
            return
        print("\"%s\" is going to be DROPPED!" % self.pathSmplHash)
        print("Print \"yes\" to continue >")
        flag = input()
        if flag != "yes":
            print("Exit: wrong input")
            return
        if self.pathSmplHash is not None:
            self.sumyCol.delete_one({"pathSmplHash": self.pathSmplHash})
        if self.fileCol is not None:
            self.fileCol.drop()
        if self.infoCol is not None:
            self.infoCol.drop()
        print("Drop Done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Upload Structure to LaspDataBase')
    parser.add_argument('-u',
                        "--user",
                        type=str,
                        default=[os.environ["USER"]],
                        nargs=1)
    parser.add_argument("-ub", "--uploadBestArc", type=str, nargs='+')
    parser.add_argument(
        "-d",
        "--download",
        action="store_true",
    )
    parser.add_argument("-a", "--auto", action="store_true")
    parser.add_argument("-k",
                        "--keys",
                        type=str,
                        nargs="*",
                        help="Add necessary description of the structure")
    args = parser.parse_args()
    print(args)


    # storage5
    # host =
    if args.uploadBestArc:
        loader = LaspDBGlobalMin("xxx", "xxx", "xxx", "xxxx")
        for file in args.uploadBestArc:
            struc = allstr()
            struc.readfile(file)
            print(file)
            loader.DBUploadStructure(
                struc[-1],
                type=0,
                author=args.user[0],
                keys=args.keys,
            )
    if args.download:
        pass
