if __name__ == "__main__":
	import pygmo as pg
	import numpy as np
	import pandas as pd
	from multi_bfe import multi_bfre, multi_bfre2
	from ExternalR import ExternalR2
	from multiprocessing.pool import ThreadPool
	#from multiprocessing import Pool
	from FireDrillUDP import FireDrillUDP3, FireDrillUDP4
	import sched, time
	from threading import Timer
	from queue import Queue
	import math

	class my_constrained_udp:
		def fitness(self, x):
			print(x)
			obj = 0
			for i in range(3):
				obj += (x[2*i-2]-3)**2 / 1000. - (x[2*i-2]-x[2*i-1]) + math.exp(20.*(x[2*i - 2]-x[2*i-1]))
			ce1 = 4*(x[0]-x[1])**2+x[1]-x[2]**2+x[2]-x[3]**2
			ce2 = 8*x[1]*(x[1]**2-x[0])-2*(1-x[1])+4*(x[1]-x[2])**2+x[0]**2+x[2]-x[3]**2+x[3]-x[4]**2
			ce3 = 8*x[2]*(x[2]**2-x[1])-2*(1-x[2])+4*(x[2]-x[3])**2+x[1]**2-x[0]+x[3]-x[4]**2+x[0]**2+x[4]-x[5]**2
			ce4 = 8*x[3]*(x[3]**2-x[2])-2*(1-x[3])+4*(x[3]-x[4])**2+x[2]**2-x[1]+x[4]-x[5]**2+x[1]**2+x[5]-x[0]
			ci1 = 8*x[4]*(x[4]**2-x[3])-2*(1-x[4])+4*(x[4]-x[5])**2+x[3]**2-x[2]+x[5]+x[2]**2-x[1]
			ci2 = -(8*x[5] * (x[5]**2-x[4])-2*(1-x[5]) +x[4]**2-x[3]+x[3]**2 - x[4])
			return [obj, ce1,ce2,ce3,ce4,ci1,ci2]
		def get_bounds(self):
			return ([-5]*6,[5]*6)
		def get_nic(self):
			return 2
		def get_nec(self):
			return 4
		def gradient(self, x):
			return pg.estimate_gradient(lambda x: self.fitness(x), x)


	MAX_NUM_PARALLEL = 12
	q = Queue(maxsize=MAX_NUM_PARALLEL)
	qwait = Queue()
	bpool = ThreadPool(processes=16)
	def batch_fitness_callback(prob, dvs):
		nx = prob.dim
		ny = 1 + prob.get_nic()
		ndvs = int(len(dvs) / nx)
		dvs_arr = np.reshape(dvs, (ndvs, nx))
		fit_arr = [0] * ndvs * ny

		def getProc(dv):
			return ExternalR2(
				cmd = "Rscript",
				args=["./optim/FireDrillUDP.R"],
				cwd="./").start_fitness(dv)


		# fill all queues:
		for idx, dv in enumerate(dvs_arr):
			if not q.full():
				er2Prom = bpool.apply_async(lambda: getProc(dv))
				q.put(item=(idx, er2Prom))
			else:
				qwait.put(item=(idx, dv))
		

		while not q.empty():
			idx, er2Prom = q.get()
			er2Prom.wait()
			er2 = er2Prom.get()
			res = er2.wait()

			# Now put the result in fit_arr:
			for yidx, y in enumerate(res):
				fit_arr[idx * ny + yidx] = y

			if not qwait.empty():
				new_idx, new_dv = qwait.get()
				new_er2Prom = bpool.apply_async(lambda: getProc(new_dv))
				q.put(item=(new_idx, new_er2Prom))
		
		return fit_arr


	udp = FireDrillUDP4()
	#udp = my_constrained_udp()
	#algo = pg.algorithm(pg.nlopt("slsqp"))
	#algo = pg.scipy_optimize(method="SLSQP", tol=1e-8)
	#algo.set_verbosity(1)

	class temp:
		def __init__(self):
			self.champs_x = []
			self.champs_f = []
			self.numLeft = 128

	t = temp()

	def addToPool(a='slsqp'):
		prob = pg.problem(udp)		
		prob.c_tol = [1e-6] * (prob.get_nic() + prob.get_nec())

		algo = None
		if a == "auglag":
			algo = pg.algorithm(uda = pg.nlopt('auglag'))
			algo.extract(pg.nlopt).local_optimizer = pg.nlopt('var2')
		elif a == "cobyla":
			#algo = pg.scipy_optimize(method='COBYLA')
			algo = pg.algorithm(uda = pg.nlopt('cobyla'))
		elif a == "slsqp":
			#algo = pg.scipy_optimize(method='SLSQP')
			algo = pg.algorithm(pg.nlopt("slsqp"))
		elif a == "krylov":
			#prob = pg.unconstrain(prob=prob, method="kuri")
			algo = pg.scipy_optimize(method='trust-krylov')
		else:
			prob = pg.unconstrain(prob=prob, method="kuri")
			algo = pg.scipy_optimize(method="BFGS")
		algo.set_verbosity(10000)

		def temp():
			pop = pg.population(prob=prob, size = 1)
			pop1 = algo.evolve(pop)
			if pop1.champion_f[0] < 5:
				print("{} terminated with champions and evals f/g:".format(a))
				print(pop1.champion_x)
				print(pop1.champion_f)
				print(pop1)
				t.champs_x.append(pop1.champion_x)
				t.champs_f.append(pop1.champion_f)
			if t.numLeft > 0:
				addToPool(a)
				t.numLeft = t.numLeft - 1

		bpool.apply_async(func=temp)
	
	for i in range(3):
		addToPool('auglag')
	for i in range(4):
		addToPool('cobyla')
	for i in range(4):
		addToPool('slsqp')
	for i in range(3):
		addToPool('krylov')
	for i in range(4):
		addToPool()
	

	def checkDone():
		while t.numLeft > 0:
			time.sleep(30)
	bpool.apply(func=checkDone)

	# isl = pg.island(
	# 	udi=pg.mp_island(),
	# 	#algo=pg.nlopt(solver="slsqp"),
	# 	algo=pg.scipy_optimize(method="SLSQP", tol=1e-8),
	# 	prob=prob,
	# 	b=pg.bfe(udbfe=multi_bfre2()),
	# 	size=MAX_NUM_PARALLEL
	# )

	
	# bpool = ThreadPool(processes=3)
	# def logger():
	# 	while True:
	# 		print(pop.champion_x)
	# 		print(pop.champion_f)
	# 		time.sleep(5)

	# bpool.apply_async(logger)



	# isl.evolve(n=10000)
	# isl.wait()

	jx = pd.DataFrame(champs_x).to_json(orient='values')
	jf = pd.DataFrame(champs_f).to_json(orient='values')

	with open("jx.json", "w") as jxf:
		jxf.write(jx)
	with open("jf.json", "w") as jff:
		jff.write(jf)

	print(5)
