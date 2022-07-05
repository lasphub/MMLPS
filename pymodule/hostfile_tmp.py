import os
import time
import multiprocessing
# from multiprocessing import Pool
import glob
import traceback
import subprocess


class Hostfile(object):
    def __init__(self, workdir, cpuperjob, setmasternode=0):
        self.workdir = workdir
        self.cpuperjob = cpuperjob
        self.setmasternode = setmasternode

    def setHostfile(self, inputhostfile=False):
        if not inputhostfile:
            hostDict = self.getProc()
        else:
            hostDict = self.getProcfromfile(inputhostfile)

        divHost, totalproc = self.alloProc(hostDict, self.cpuperjob)
        print(divHost)
        # get the pool size
        poolsize = len(divHost)
        print(poolsize)
        # dump the hostfiles
        self.hostInfo = self.dumpHost(divHost)
        return self.hostInfo, poolsize, totalproc

    def getProcfromfile(self, hostfile):

        hostdict = {}
        if self.setmasternode != 0:
            first = True
            print('set masternode')
        else:
            first = False

        # read hostfile
        fp = open(hostfile, "r")
        for x in fp:
            if first:
                masternode = x.split()[0]
                first = False
                continue
            line = x.split()
            hostdict[line[0]] = int(line[1])
        fp.close()

        return hostdict

    def getProc(self):
        hostdict = {}
        if self.setmasternode != 0:
            first = True
            print('set masternode')
        else:
            first = False

        # get environment variable of host file list in group cluster
        hostfile = os.environ["PE_HOSTFILE"]
        print(hostfile)
        fp = open(hostfile, "r")
        for x in fp:
            if first:
                masternode = x.split()[0]
                first = False
                continue
            line = x.split()
            hostdict[line[0]] = int(line[1])
        fp.close()
        #        if self.setmasternode: hostdict.pop(masternode)

        # the following part are designed for fudannhpcc system
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

        f = open('hostfile_now', 'w')
        for key, val in hostdict.items():
            f.write('%s    %s\n' % (key, val))
        f.close()

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

        print("Availiable proc number: "), totProc
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


def runprog_cluster(workdir, prog, ncpus, hostInfo, rootdir, poolcount=0):

    try:
        os.chdir(rootdir)
        cwd = os.getcwd()
        exit = False
        # print multiprocessing.current_process().name
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
        # totalrun = 'source ~/.bashrc; '+ mpiprog
        totalrun = mpiprog
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
        fout.close()
        return
    except Exception as e:
        traceback.print_exc()
        print('error')
        raise e


def runprog_cluster_manual(workdir,
                           prog,
                           hostInfo,
                           ncpus,
                           rootdir,
                           env,
                           poolcount=0,
                           maxtime=False,
                           Lcheck=False):

    try:
        os.chdir(rootdir)
        cwd = os.getcwd()
        exit = False
        # print multiprocessing.current_process().name
        nodeInfo = hostInfo[
            int(multiprocessing.current_process().name.split("-")[-1]) - 1 -
            poolcount]
        os.chdir(workdir)
        mf = os.path.join(cwd, nodeInfo[0])
        # mpiprog = "/opt/intel/impi/5.0.2.044/intel64/bin/mpirun -machinefile %s -np %d "%(mf, ncpus) + prog
        # mpiprog = "mpirun -machinefile %s -env I_MPI_DEVICE rdma:OpenIB-cma -np %d "%(mf, ncpus) + prog
        #        mpiprog = "source /home2/shang/.bashrc; mpirun -machinefile %s -np %d "%(mf, ncpus) + prog
        mpiprog = " mpirun  -rsh=ssh -machinefile %s -np %d " % (mf,
                                                                 ncpus) + prog
        #       mpiprog = "mpirun -machinefile %s -np %d "%(mf, ncpus) + prog
        fout1 = open("proginfo", "w")
        fout1.write("Current process: " +
                    multiprocessing.current_process().name + "\n")
        fout1.write("Host file: " + mf + "\n")
        fout = open('output', 'w')
        # totalrun = 'source ~/.bashrc; '+ mpiprog
        totalrun = mpiprog
        # nodename = os.popen('head -1 %s'%mf).readline().strip()
        # totalrun = 'source ~/.bashrc; ssh '+nodename +"; "+mpiprog
        fout1.write(totalrun + "\n")
        child = subprocess.Popen(mpiprog,
                                 stdout=fout,
                                 stderr=fout,
                                 shell=True,
                                 executable='/bin/bash',
                                 preexec_fn=os.setpgrp)
        # print ('start run job in %s'%workdir)
        fout1.write('pid   %d\n' % child.pid)
        fout1.close()
        pid = child.pid
        if maxtime:
            alltime = 0
        while not exit:
            time.sleep(10)
            returnCode = child.poll()
            if glob.glob('killsignal'):
                os.kill(-pid, 9)
                time.sleep(3)
                fout.write('kill %s\n' % pid)
                a = os.waitpid(pid, 0)
                print(a)
                time.sleep(3)
                exit = True

            if Lcheck:
                ff = open('output', 'r')
                for Line in ff:
                    if 'On entry to' in Line:
                        os.kill(-pid, 9)
                        time.sleep(3)
                        fout.write('kill %s for problem converge\n' % pid)
                        a = os.waitpid(pid, 0)
                        print(a)
                        time.sleep(3)
                        exit = True
                        break
                ff.close()

            if isinstance(returnCode, int):
                if returnCode == 0:
                    fout.write('successfully done\n')
                else:
                    fout.write('something wrong: returnCode  %d\n' %
                               returnCode)
                exit = True
            if maxtime:
                alltime = alltime + 10
                if alltime > maxtime:
                    os.kill(-pid, 9)
                    os.system('pkill -9 %s' % prog)
                    time.sleep(3)
                    fout.write('time out :kill %s\n' % pid)
                    a = os.waitpid(pid, 0)
                    print(a)
                    time.sleep(3)
                    exit = True

        # fout1.write(str(child.poll())+'\n')
        fout.write('exit\n')
        fout.close()
        return
    except Exception as e:
        traceback.print_exc()
        raise e


if __name__ == "__main__":

    ncpu = 24
    rootdir = os.getcwd()
    Host = Hostfile(rootdir, ncpu, 0)
    hostInfo, poolsize, totalproc = Host.setHostfile()
