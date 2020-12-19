import pygmo as pg



if __name__ == "__main__":

    class sphere_1d:
        def get_bounds(self):
            return ([0], [1])
        
        def fitness(self, dv):
            i = 0
            while (i < 1e7): i = i +1

            return [dv[0]**2]



    mpbfe = pg.mp_bfe(chunksize=1)
    pg.mp_bfe.init_pool(processes=4)
    bfe = pg.bfe(udbfe=mpbfe)

    udp = pg.schwefel(dim = 19)
    prob = pg.problem(udp)
    algo = pg.algorithm(pg.simulated_annealing())
    pop = pg.population(prob=prob, b=bfe, size=22)

    archi = pg.archipelago(n=1, algo=algo, pop=pop)

    print(archi.get_champions_f())
    archi.evolve()
    archi.wait()
    print(archi.get_champions_f())

    print(pop.champion_f)
    pop = algo.evolve(pop)
    print(pop.champion_f)

    #ils = pg.island(algo=algo, udi=pg.thread_island(), pop=pop)

    print(5)

    #temp = ils.evolve(n=10)
    #ils.wait()

    print(7)
