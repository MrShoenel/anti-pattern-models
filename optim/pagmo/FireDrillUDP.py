import pygmo as pg
from ExternalR import ExternalR, ExternalR2
#from queue import Queue



class FireDrillUDP:

	def __init__(self, callback, bfeCallback, dim=23, nic=10):
		self.dim = dim
		self.nic = nic
		self.callback = callback
		self.bfeCallback = bfeCallback

	def fitness(self, x):
		return self.callback(x)
	
	def batch_fitness(self, dvs):
		return self.bfeCallback(dvs)

#	def gradient(self, x):
#		return pg.estimate_gradient(lambda x: self.fitness(x), x) 

	def get_nec(self):
		return 0

	def get_nic(self):
		return self.nic

	def get_bounds(self):
		return ([0] * self.dim, [1] * self.dim)


class FireDrillUDP2:
	def __init__(self, dim=23, nic=10):
		self.dim = dim
		self.nic = nic
	
	def gradient(self, x):
		return pg.estimate_gradient(lambda x: self.fitness(x), x) 

	def feasibility_x(self, x):
		b1 = x[0]
		b2 = x[1]
		b3 = x[2]

		if b1 < .01 or b3 > .99:
			return False
		if ((b2 - b1) < .025) or ((b3 - b2) < .025):
			return False
		return True

	def fitness(self, x):
		e = ExternalR(
			cmd = "Rscript",
			args=["./optim/FireDrillUDP.R"],
			cwd="./")
		res = e.fitness(x=x)
		#print(res[0])
		#print((x, res))
		# if res[0] < 1e308:
		# 	print(res[0])
		return res

	def get_nec(self):
		return 0

	def get_nic(self):
		return self.nic

	def get_bounds(self):
		return ([0] * self.dim, [1] * self.dim)



class FireDrillUDP3:
	def __init__(self, callback, dim=23, nic=10):
		self.dim = dim
		self.nic = nic
		self.callback = callback
	
	def gradient(self, x):
		return pg.estimate_gradient(lambda x: self.fitness(x), x) 

	def feasibility_x(self, x):
		b1 = x[0]
		b2 = x[1]
		b3 = x[2]

		if b1 < .01 or b3 > .99:
			return False
		if ((b2 - b1) < .025) or ((b3 - b2) < .025):
			return False
		return True

	def fitness(self, x):
		return self.batch_fitness(dvs=x)
	
	def batch_fitness(self, dvs):
		return self.callback(self, dvs)
		# nx = self.get_nx()
		# ny = 1 + self.nic
		# ndvs = int(len(dvs) / nx)
		# dvs_arr = np.reshape(dvs, ndvs, nx)
		# fit_arr = [0] * ndvs * ny

		# def getProc():
		# 	return ExternalR2(
		# 		cmd = "C:\\Program Files\\R\\R-4.0.2\\bin\\Rscript.exe",
		# 		args=["FireDrillUDP.R"],
		# 		cwd="C:\\repos\\lnu_anti-patterns\\optim")


		# # fill all queues:
		# for idx, dv in enumerate(dvs_arr):
		# 	if not self.q.full():
		# 		er2 = getProc()
		# 		self.q.put(item=(idx, er2.start_fitness(dv)))
		# 	else:
		# 		self.qwait.put(item=(idx, dv))
		

		# while not self.q.empty():
		# 	idx, er2 = self.q.get()
		# 	res = er2.wait()

		# 	# Now put the result in fit_arr:
		# 	for yidx, y in enumerate(res):
		# 		fit_arr[idx * ny + yidx] = y

		# 	if not self.qwait.empty():
		# 		new_idx, new_dv = self.qwait.get()
		# 		new_er2 = getProc()
		# 		self.q.put(item=(new_idx, new_er2.start(new_dv)))
		
		# return fit_arr


	def get_nec(self):
		return 0

	def get_nic(self):
		return self.nic

	def get_bounds(self):
		return ([0] * self.dim, [1] * self.dim)




class FireDrillUDP4:
	def __init__(self, dim=23, nic=10):
		self.dim = dim
		self.nic = nic
	
	def gradient(self, x):
		return pg.estimate_gradient(lambda x: self.fitness(x), x) 

	def feasibility_x(self, x):
		b1 = x[0]
		b2 = x[1]
		b3 = x[2]

		if b1 < .01 or b3 > .99:
			return False
		if ((b2 - b1) < .025) or ((b3 - b2) < .025):
			return False
		return True

	def fitness(self, x):
		# print(x)
		e = ExternalR(
			cmd = "Rscript",
			args=["./optim/FireDrillUDP.R"],
			cwd="./")
		res = e.fitness(x=x)
		#print(res)
		return res

	def get_nec(self):
		return 0

	def get_nic(self):
		return self.nic

	def get_bounds(self):
		return ([0] * self.dim, [1] * self.dim)
