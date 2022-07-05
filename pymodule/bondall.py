from functools import reduce


class bondtype(object):
    def __init__(self, bondnum):
        self.num = bondnum
        self.order = bondnum % 10


class bondpredict(object):
    def __init__(self, allbondlist):
        self.allbond = []
        for bond in allbondlist:
            self.allbond.append(bondtype(bond))
        self.allbond.sort(key=lambda x: x.num, reverse=True)
        self.sortbyorder()

        self.allposible = []
        # self.setorder()

    def sortbyorder(self):
        self.allorder = []
        for bond in self.allbond:
            if bond.order not in self.allorder:
                self.allorder.append(bond.order)
        self.allorder.sort(reverse=True)
        print(self.allorder)

    def reduce_bondcombine2(self, nowbond, numneed):
        orderneed = numneed - sum(nowbond)

        for order in self.allorder:
            if (order <= nowbond[-1] and order <= orderneed):

                temp1 = nowbond[:]
                temp1.append(order)
                orderneed = orderneed - order

                if nowbond[-1] != 1:
                    if orderneed == 0:
                        if temp1[-1] == 1:
                            self.allposible.append(tuple(temp1))
                            nowbond[-1] = nowbond[-1] - 1
                            return self.reduce_bondcombine(nowbond, numneed)
                        else:
                            self.allposible.append(tuple(temp1))
                            temp2 = temp1[:]
                            temp2[-1] = temp2[-1] - 1
                            nowbond[-1] = nowbond[-1] - 1
                            return self.reduce_bondcombine(
                                nowbond, numneed), self.reduce_bondcombine(
                                    temp2, numneed)

                    elif orderneed != 0:
                        if temp1[-1] == 1:
                            nowbond[-1] = nowbond[-1] - 1
                            return self.reduce_bondcombine(
                                temp1, numneed), self.reduce_bondcombine(
                                    nowbond, numneed)
                        else:
                            temp2 = temp1[:]
                            temp2[-1] = temp2[-1] - 1
                            nowbond[-1] = nowbond[-1] - 1
                            return self.reduce_bondcombine(
                                temp1, numneed), self.reduce_bondcombine(
                                    nowbond, numneed
                                )  # ,self.reduce_bondcombine(temp2,numneed)

                else:
                    if orderneed == 0:
                        self.allposible.append(tuple(temp1))
                        return
                    elif orderneed != 0:
                        return self.reduce_bondcombine(temp1, numneed)

    def reduce_bondcombine(self, nowbond, numneed):
        orderneed = numneed - sum(nowbond)
        if orderneed == 0:  # for onebond situation
            self.allposible.append(tuple(nowbond))
            nowbond[-1] = nowbond[-1] - 1
            return self.reduce_bondcombine(nowbond, numneed)

        for order in self.allorder:
            if (order <= nowbond[-1] and order <= orderneed):

                temp1 = nowbond[:]
                temp1.append(order)
                orderneed = orderneed - order

                if nowbond[-1] != 1:
                    if orderneed == 0:
                        if temp1[-1] == 1:
                            self.allposible.append(tuple(temp1))
                            nowbond[-1] = nowbond[-1] - 1
                            return self.reduce_bondcombine(nowbond, numneed)
                        else:
                            self.allposible.append(tuple(temp1))
                            temp2 = temp1[:]
                            temp2[-1] = temp2[-1] - 1
                            nowbond[-1] = nowbond[-1] - 1
                            return self.reduce_bondcombine(
                                nowbond, numneed), self.reduce_bondcombine(
                                    temp2, numneed)

                    elif orderneed != 0:
                        if temp1[-1] == 1:  # impossible situation
                            nowbond[-1] = nowbond[-1] - 1
                            print('error?')
                            return self.reduce_bondcombine(
                                temp1, numneed), self.reduce_bondcombine(
                                    nowbond, numneed)
                        else:
                            # temp2= temp1[:]
                            # temp2[-1]= temp2[-1]-1
                            nowbond[-1] = nowbond[-1] - 1
                            return self.reduce_bondcombine(
                                temp1, numneed), self.reduce_bondcombine(
                                    nowbond, numneed
                                )  # ,self.reduce_bondcombine(temp2,numneed)

                else:
                    if orderneed == 0:
                        self.allposible.append(tuple(temp1))
                        return
                    elif orderneed != 0:
                        return self.reduce_bondcombine(temp1, numneed)

#                    if orderneed == 0 and order == 1:
#                        return temp1
#                    elif (orderneed ==0 and order != 1):
#                        ordernew = order -1
#                        orderneednew = orderneed +1
#                        temp2.append(ordernew)
#                        return temp1,self.reduce_bondcombine(temp2,orderneednew)
#                    elif orderneed !=0 and order != 1:
#                        ordernew = order -1
#                        orderneednew = orderneed +1
#                        temp2.append(ordernew)
#                        return self.reduce_bondcombine(temp1,orderneed),self.reduce_bondcombine(temp2,orderneednew)
#                    elif orderneed != 0 and order == 1:
#                        return self.reduce_bondcombine(temp1,orderneed)

    def grepbondbyorder(self, order):
        templist = []
        for bond in self.allbond:
            if bond.order == order:
                templist.append(bond)
        return templist

    def fillbondtype(self, order, nbond, tempin, templist):
        ibond = len(templist)

        for i, bond in enumerate(templist):
            if bond.num == tempin[-1].num:

                if len(tempin) == nbond:  # for nbond ==1
                    if (i == ibond - 1):
                        self.orderfill.append([x.num for x in tempin])
                        return
                    if (i < ibond - 1):
                        self.orderfill.append([x.num for x in tempin])
                        tempin[-1] = templist[i + 1]
                        return self.fillbondtype(order, nbond, tempin,
                                                 templist)

                tempfill = tempin[:]
                tempfill.append(bond)
                if len(tempfill) == nbond and i < ibond - 1:
                    self.orderfill.append([x.num for x in tempfill])
                    tempfill2 = tempin[:]
                    tempfill2.append(templist[i + 1])
                    tempin[-1] = templist[i + 1]
                    return self.fillbondtype(order, nbond, tempfill2,
                                             templist), self.fillbondtype(
                                                 order, nbond, tempin,
                                                 templist)
                elif len(tempfill) == nbond and i == (ibond - 1):
                    self.orderfill.append([x.num for x in tempfill])
                    return
                elif len(tempfill) < nbond and i < (ibond - 1):
                    tempin[-1] = templist[i + 1]
                    return self.fillbondtype(order, nbond, tempin,
                                             templist), self.fillbondtype(
                                                 order, nbond, tempfill,
                                                 templist)
                elif len(tempfill) < nbond and i == (ibond - 1):
                    return self.fillbondtype(order, nbond, tempfill, templist)

    def filling(self, orderlist):
        self.bondcombine = []

        temporder = orderlist[0]
        nbond = 0
        for i in range(len(orderlist)):
            if (temporder != orderlist[i]):
                self.orderfill = []
                tempin = []
                templist = self.grepbondbyorder(temporder)
                tempin.append(templist[0])
                self.fillbondtype(temporder, nbond, tempin, templist)
                self.bondcombine.append(self.orderfill)
                nbond = 1
                temporder = orderlist[i]
            else:
                nbond = nbond + 1

        self.orderfill = []
        tempin = []
        templist = self.grepbondbyorder(temporder)
        tempin.append(templist[0])
        self.fillbondtype(temporder, nbond, tempin, templist)
        self.bondcombine.append(self.orderfill)

    def combinefrag(self, list1, list2):
        list3 = []
        i = 0
        for frag1 in list1:
            for frag2 in list2:
                list3.append([])
                list3[i] = frag1 + frag2
                i = i + 1
        return list3

    def bondcombine(self, num):
        test = []
        test.append(self.allorder[0])
        self.reduce_bondcombine(test, num)
        # for orderlist in self.allposible:
        # self.filling(orderlist)

        finalresult = []
        for orderlist in self.allposible:
            self.filling(orderlist)
            # fragnum = len(self.bondcombine)
            allcombine = reduce(self.combinefrag, self.bondcombine)
            finalresult.extend(allcombine)

        return finalresult


if __name__ == "__main__":
    bond = bondpredict([13, 103, 12, 102, 11, 101, 1])
    allbond = bond.bondcombine(4)
    print(allbond)
