import pygmo as pg
from multi_bfe import multi_bfre
from ExternalR import ExternalR
from multiprocessing.pool import ThreadPool
#from multiprocessing import Pool
from FireDrillUDP import FireDrillUDP
import sched, time
from threading import Timer

tpool = ThreadPool(processes=20)
bpool = ThreadPool(processes=20)


def bfeCallback(x):
    def temp():
        return udp.fitness(x)
    return bpool.apply_async(func=temp)

muBfe = multi_bfre(callback=bfeCallback)    
bfe = pg.bfe(udbfe = muBfe)

def fdCallback(x):
    def temp():
        e = ExternalR(
            cmd = "C:\\Program Files\\R\\R-4.0.2\\bin\\Rscript.exe",
            args=["FireDrillUDP.R"],
            cwd="C:\\repos\\lnu_anti-patterns\\optim")
        return e.fitness(x=x)
    res = tpool.apply(func=temp)
    return res

def fdBfeCallback(dvs):
    muBfe(prob=prob, dvs=dvs)

udp = FireDrillUDP(callback=fdCallback, bfeCallback=fdBfeCallback)
prob = pg.problem(udp)
algo = pg.scipy_optimize(method="SLSQP", tol=1e-8)
pop = pg.population(prob, 16, b=bfe)
algo.evolve(pop)

print(pop.champion_f)
print(pop.problem.get_fevals())
print(pop.problem.get_gevals())
print(pop.champion_x)

print(5)
