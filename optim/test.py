import pygmo as pg
from pygmo import problem, algorithm, population, de, de1220



class udp_mlmrc:
	def __init__(self):
		pass

	def fitness(self, x):
		# Here, we pass x to the R process
		pass


class add_gradient:
    def __init__(self, prob):
            self.prob = pg.problem(prob)
    def fitness(self, x):
        return self.prob.fitness(x)
    def get_bounds(self):
        return self.prob.get_bounds()
    def get_nec(self):
        return self.prob.get_nec()
    def get_nic(self):
        return self.prob.get_nic()
    def get_nobj(self):
        return self.prob.get_nobj()
    def gradient(self, x):
        return pg.estimate_gradient(lambda x: self.fitness(x), x)



class sphere_1d:
	def get_bounds(self):
		return ([0], [1])
	
	def fitness(self, dv):
		return [dv[0]**2]


prob = problem(add_gradient(sphere_1d()))
algo = pg.algorithm(uda = pg.mbh(pg.nlopt("slsqp"), 20, .2))
pop = pg.population(prob, 1)
# Solve the problem
new_pop = algo.evolve(pop) 
# Collect information
print(new_pop.champion_f) 
print(new_pop.problem.get_fevals())
print(new_pop.problem.get_gevals())

print("---------")

#prob = sphere_1d()
#algo = algorithm(de(500))
algo = algorithm(de1220(gen = 500))
pop = population(prob, 20)
new_pop = algo.evolve(pop)
print(new_pop.champion_f)
print(new_pop.problem.get_fevals())
print(new_pop.problem.get_gevals())
