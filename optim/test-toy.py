import time
import random


def f_log_decor(orig_fitness_function):
    def new_fitness_function(self, dv):
        if hasattr(self, "dv_log"):
            self.dv_log.append(dv)
        else:
            self.dv_log = [dv]
        return orig_fitness_function(self, dv)
    return new_fitness_function


class toy_problem:
    def __init__(self, dim):
        self.cnt = 0
        self.dim = dim

    def fitness(self, x):
        return [sum(x), 1 - sum(x*x), - sum(x) + .1]

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



if __name__ == "__main__":
    import pygmo as pg
    a_cstrs_sa = pg.algorithm(pg.cstrs_self_adaptive(iters=1000))
    #a_cstrs_sa = pg.algorithm(pg.de1220(gen = 500))
    #a_cstrs_sa = pg.cstrs_self_adaptive(algo=pg.de1220(gen = 500))
    tp = toy_problem(50)
    p_toy = pg.problem(tp)
    #p_toy = pg.problem(pg.decorator_problem(tp, fitness_decorator=f_log_decor))
    p_toy.c_tol = [1e-7, 1e-7]
    archi = pg.archipelago(n=32,algo=a_cstrs_sa, prob=p_toy, pop_size=70)
    print(archi) 

    print(archi.get_champions_f())

    archi.evolve() 
    print("------")
    print(tp.cnt)
    print(p_toy.get_fevals())
    print(p_toy.get_gevals())
    print(p_toy.get_hevals())
    print(archi) 


    archi.wait()
    print("------")
    print(tp.cnt)
    print(p_toy.get_fevals())
    print(p_toy.get_gevals())
    print(p_toy.get_hevals())
    archi.wait_check()
    print(archi.get_champions_f())
