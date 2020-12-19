import numpy as np

"""This class is a user-defined batch fitness evaluator (UDBFE)
that allows mapping one or more decision vectors to singular
calls to fitness(x) and returns their results.
"""
class multi_bfre:
	def __init__(self, callback):
		self.callback = callback

	def __call__(self, prob, dvs):
		nx = prob.get_nx()
		dvs_arr = np.reshape(dvs, (int(len(dvs) / nx), nx))
		promises = [self.callback(x) for x in dvs_arr]
		for p in promises:
			p.wait()
		
		res = [p.get() for p in promises]
		# Remember that fitness() returns a vector
		return np.array(res).flatten()



class multi_bfre2:
	def __call__(self, prob, dvs):
		return prob.batch_fitness(dvs)
