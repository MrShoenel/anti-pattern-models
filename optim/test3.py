import pygmo as pg


if __name__ == "__main__":
    #prob = pg.problem(pg.rosenbrock(dim = 30))
    udp = pg.schwefel(dim = 19)
    prob = pg.problem(udp)
    #pop1 = pg.population(prob, size=73)

    #algo = pg.algorithm(pg.sade(gen=500))
    #algo = pg.algorithm(pg.scipy_optimize(method="Nelder-Mead"))
    algo = pg.algorithm(pg.simulated_annealing())
    #algo.set_verbosity(10)
    for i in range(3):
        algo.set_verbosity(1)
        pop = pg.population(prob=prob, size=22)
        pop = algo.evolve(pop)
        print(pop.champion_f)
    

    
    #archi = pg.archipelago(n=4,algo=algo, pop=pop1)
    print(archi) 
    archi.evolve() 
    archi.wait()
    archi.wait_check()
    print(archi)




import pygmo as pg
# The user-defined problem
udp = pg.schwefel(dim = 20)
# The pygmo problem
prob = pg.problem(udp)

# For a number of generation based algorithms we can use a similar script to run and average over 25 runs.
udas = [pg.sade(gen=500), pg.de(gen=500), pg.de1220(gen=500), pg.pso(gen=500), pg.bee_colony(gen=250, limit=20)]
for uda in udas: 
    logs = []
    for i in range(25):
        algo = pg.algorithm(uda)
        algo.set_verbosity(1) # regulates both screen and log verbosity
        pop = pg.population(prob, 20)
        pop = algo.evolve(pop)
        logs.append(algo.extract(type(uda)).get_log())
    logs = np.array(logs)
    avg_log = np.average(logs,0)
    plt.plot(avg_log[:,1],avg_log[:,2]-418.9829*20 , label=algo.get_name())
