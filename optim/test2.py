import pygmo as pg
import numpy as np
from pygmo import problem, algorithm, population, de, de1220


class my_udp:
    def fitness(self, x):
        return (np.sin(x[0]+x[1]-x[2]), x[0] + np.cos(x[2]*x[1]), x[2])
    def get_bounds(self):
        return ([-1,-1,-1],[1,1,1])
    def get_nec(self):
        return 1
    def get_nic(self):
        return 1


prob = problem(my_udp())
pop = pg.population(prob, 16, seed=1337)
#algo = pg.algorithm(pg.scipy_optimize(method="SLSQP"))
algo = pg.algorithm(pg.de1220())
algo.set_verbosity(10)
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
