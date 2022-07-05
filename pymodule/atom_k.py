import PeriodicTable as PT


class S_atom(object):
    def __init__(self, coord, element):
        self.xyz = coord
        self.ele = element
        self.id = 0
        self.elesymbol = PT.Eletable[self.ele - 1]
        self.force = []
        self.symf = []
        self.species = 0
        self.bondlist = []
        self.bondtype = []
        self.Ctye = 'empty'
        self.charge = 0
