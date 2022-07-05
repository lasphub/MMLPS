#! /home10/shiyf/bin/miniconda3/bin/python
import os


class LaspInDict(dict):
    def __init__(self, filename):
        super(dict, self).__init__()

        k = ""
        v = 0
        blockFlag = False

        with open(filename, "r") as fp:
            lines = fp.readlines()
        i = -1
        while i < len(lines) - 1:
            i += 1
            li = lines[i].strip()
            if li == "" or "#" in li:
                continue
            if blockFlag is True:
                if li.startswith("%block"):
                    blockFlag = False
                    self[k.lower()] = v
                else:
                    v.extend([float(x) for x in li.split()])
            else:
                if li.startswith("%block"):
                    blockFlag = True
                    k = li.split()[1]
                    v = []
                else:
                    k = str(li.split()[0])
                    try:
                        v = float(li.split()[1])
                    except ValueError:
                        v = str(li.split()[1])
                    self[k.lower()] = v


globalDict = {}


def ReadLaspIn():
    """
    Read Global Variables from 'input' or 'lasp.in' file.
    All variables will be stored in globalDict while all keys are in lowwercase
    """
    fname = None
    if os.path.exists("./input"):
        fname = "input"
    elif os.path.exists("./lasp.in"):
        fname = "lasp.in"
    if fname:
        globals()["globalDict"] = LaspInDict(fname)
    else:
        print("Please write lasp.in or input files!")


if __name__ == "__main__":
    ReadLaspIn()
    print(globalDict)
