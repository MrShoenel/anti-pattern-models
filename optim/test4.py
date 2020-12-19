import pygmo as pg
import numpy as np
import subprocess as sb
import multiprocessing as mp
from multiprocessing.pool import ThreadPool
from joblib import Parallel, delayed



# tpool = ThreadPool(processes=3)

# def testFunc():
#	 i = 0
#	 while (i < 1e7):
#		 i = i+1
#	 return i

# res1 = tpool.apply_async(func=testFunc)
# res2 = tpool.apply_async(func=testFunc)
# res1.wait()
# res2.wait()

# sb.call(args=[], executable='')


class ext_proc_objf:
	def __init__(self, pathToExe):
		self.pathToExe = pathToExe



class my_bfre:
	def __call__(self, prob, dvs):
		nx = prob.get_nx()
		dvs_arr = np.reshape(dvs, (int(len(dvs) / nx), nx))
		res = [prob.fitness(x) for x in dvs_arr]
		# Remember that fitness() returns a vector
		return np.array(res).flatten()


bfe = pg.bfe(udbfe = my_bfre())


class sphere_1d:
	def get_bounds(self):
		return ([0], [1])
	
	def fitness(self, dv):
		return [dv[0]**2]

	def gradient(self, x):
		return pg.estimate_gradient(lambda x: self.fitness(x), x)


class toy_problem:
	def __init__(self, dim):
		self.cnt = 0
		self.dim = dim

	def fitness(self, x):
        # All equality constraints are in the form g(x)=0,
        # while inequalities are in the form g(x)<=0 
		return [sum(x), sum(x*x) - 1, -1*(sum(x) - .1)]

	def gradient(self, x):
		return pg.estimate_gradient(lambda x: self.fitness(x), x) 

	def get_nec(self):
		return 1

	def get_nic(self):
		return 1

	def get_bounds(self):
		return ([-1] * self.dim, [1] * self.dim)

	def get_name(self):
		return "A toy problem"

	def get_extra_info(self):
		return "\tDimensions: " + str(self.dim)




#udp = sphere_1d()
udp = toy_problem(dim=159)
prob = pg.problem(udp)
#algo = pg.algorithm(uda = pg.mbh(pg.nlopt("slsqp"), 19, .2))
#algo = pg.algorithm(pg.cstrs_self_adaptive(algo=pg.nlopt("slsqp")))
#algo = pg.algorithm(pg.cstrs_self_adaptive(iters=1000))
# WORKS
algo = pg.scipy_optimize(method="SLSQP", tol=1e-10)
#algo = pg.algorithm(pg.cstrs_self_adaptive(iters=1000))
#algo = pg.algorithm(pg.cstrs_self_adaptive(algo=pg.scipy_optimize(method="SLSQP")))

pop = pg.population(prob, 20, b=bfe)
new_pop = algo.evolve(pop)
print(new_pop.champion_f) 
print(new_pop.problem.get_fevals())
print(new_pop.problem.get_gevals())
print(new_pop.champion_x) 
