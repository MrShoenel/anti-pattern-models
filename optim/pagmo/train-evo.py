if __name__ == "__main__":
	import pygmo as pg
	import pandas as pd
	import numpy as np
	from multi_bfe import multi_bfre2
	from ExternalR import ExternalR2
	from multiprocessing.pool import ThreadPool
	#from multiprocessing import Pool
	from FireDrillUDP import FireDrillUDP, FireDrillUDP3
	import sched, time
	from threading import Timer
	from queue import Queue



	MAX_NUM_PARALLEL = 15
	q = Queue(maxsize=MAX_NUM_PARALLEL)
	qwait = Queue()
	def batch_fitness_callback(prob, dvs):
		nx = prob.dim
		ny = 1 + prob.get_nic()
		ndvs = int(len(dvs) / nx)
		dvs_arr = np.reshape(dvs, (ndvs, nx))
		fit_arr = [0] * ndvs * ny

		def getProc():
			return ExternalR2(
				cmd = "Rscript",
				args=["./optim/FireDrillUDP.R"],
				cwd="./")


		# fill all queues:
		for idx, dv in enumerate(dvs_arr):
			if not q.full():
				er2 = getProc()
				q.put(item=(idx, er2.start_fitness(dv)))
			else:
				qwait.put(item=(idx, dv))
		

		while not q.empty():
			idx, er2 = q.get()
			res = er2.wait()

			# Now put the result in fit_arr:
			for yidx, y in enumerate(res):
				fit_arr[idx * ny + yidx] = y

			if not qwait.empty():
				new_idx, new_dv = qwait.get()
				new_er2 = getProc()
				q.put(item=(new_idx, new_er2.start_fitness(new_dv)))
		
		return fit_arr


	# bpool = ThreadPool(processes=16)
	# def bfeCallback(x):
	#	 def temp():
	#		 return udp.fitness(x)
	#	 return bpool.apply_async(func=temp)

	# muBfe = multi_bfre(callback=bfeCallback)	
	# bfe = pg.bfe(udbfe = muBfe)


	udp = FireDrillUDP3(callback=batch_fitness_callback)
	prob = pg.problem(udp)
	#prob = pg.problems.unconstrain(prob=udp, method='kuri')


	# isl1 = pg.island(
	#	 udi=pg.mp_island(),
	#	 algo=pg.algorithm(pg.cstrs_self_adaptive(iters=500, algo=pg.pso())),
	#	 prob=prob,
	#	 size=7)
	# isl2 = pg.island(
	#	 udi=pg.mp_island(),
	#	 algo = pg.algorithm(pg.cstrs_self_adaptive(iters=500, algo=pg.simulated_annealing())),
	#	 prob=prob,
	#	 size=7)




	archi = pg.archipelago(
		n=0,
		t=pg.topology(udt=pg.fully_connected(n=80)),
		# udi=pg.mp_island(),
		# algo=pg.algorithm(pg.cstrs_self_adaptive(iters=1000, algo=pg.pso())),
		# prob=prob,
		# pop_size=16,
		b=pg.bfe(udbfe=multi_bfre2())
	)


	def addIsland(arc, algo):
		isl = pg.island(
			udi=pg.mp_island(),
			algo=algo,
			prob=prob,
			#b=pg.mp_bfe(),
			b=pg.bfe(udbfe=multi_bfre2()),
			size=16)
		arc.push_back(isl)
	
	for i in range(2):
		addIsland(arc=archi, algo=pg.algorithm(pg.cstrs_self_adaptive(iters=1000, algo=pg.de(gen=50))))
		addIsland(arc=archi, algo=pg.algorithm(pg.cstrs_self_adaptive(iters=1000, algo=pg.pso(gen=50))))
		addIsland(arc=archi, algo=pg.algorithm(pg.cstrs_self_adaptive(iters=1000, algo=pg.simulated_annealing())))
		addIsland(arc=archi, algo=pg.algorithm(pg.cstrs_self_adaptive(iters=1000, algo=pg.sea(gen=50))))


	bpool = ThreadPool(processes=3)
	def logger():
		#while True:
		print(archi.get_champions_x())
		print(sorted(list(map(lambda a: a[0], archi.get_champions_f()))))
		#print(archi.get_champions_f())

		jx = pd.DataFrame(archi.get_champions_x()).to_json(orient='values')
		jf = pd.DataFrame(archi.get_champions_f()).to_json(orient='values')

		with open("jx.json", "w") as jxf:
			jxf.write(jx)
		with open("jf.json", "w") as jff:
			jff.write(jf)
		
		#time.sleep(10)

	#bpool.apply_async(logger)


	for k in range(10000):
		archi.evolve()
		archi.wait()
		print(sorted(list(map(lambda a: a[0], archi.get_champions_f()))))
		logger()
		#archi.wait_check()
	

	print(archi.get_champions_x())
	print(archi.get_champions_f())

	print(5)
