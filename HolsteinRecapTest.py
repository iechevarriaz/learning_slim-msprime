import msprime
import math
import pyslim
import time

startFloat = time.process_time()
start = str(startFloat)
start2Float = time.perf_counter()
start2 = str(start2Float)
# Times are provided in years, so we convert into generations.
generation_time = 5
st = 20553234
T_0= st / generation_time
T_1 = (st + 5) / generation_time
T_2 = (st + 20) / generation_time
T_3 = (st + 35) / generation_time
T_4 = (st + 65) / generation_time
T_5 = (st + 95) / generation_time
T_6 = (st + 125) / generation_time
T_7 = (st + 775) / generation_time
T_8 = (st + 2275) / generation_time
T_9 = (st + 3275) / generation_time
T_10 = (st + 8775) / generation_time
T_11 = (st + 11775) / generation_time
T_12 = (st + 16775) / generation_time
T_13 = (st + 165775) / generation_time
T_14 = (st + 4665775) / generation_time

N_EU0 = 90
N_A0 = 120
r_EU = 0.111111
r_AS = 0.361111
#N_EU = N_EU0 / math.exp(-r_EU * T_1)
#N_A = N_A0 / math.exp(-r_AS * T_2)

population_configurations = [
    msprime.PopulationConfiguration(
        sample_size=10, initial_size=90),
    msprime.PopulationConfiguration(
        sample_size=10, initial_size=100),

]

demographic_events = [

    msprime.PopulationParametersChange(
        time=T_1, initial_size=90, growth_rate=-0.095894024, population_id=0),
    msprime.PopulationParametersChange(
        time=T_2, growth_rate=-0.24465639, population_id=0),
    msprime.PopulationParametersChange(
        time=T_3, growth_rate=-0.0560787, population_id=0),
    msprime.PopulationParametersChange(
        time=T_4, growth_rate=-0.1749704, population_id=0),
    msprime.PopulationParametersChange(
        time=T_5, growth_rate=-0.0675775, population_id=0),
    msprime.PopulationParametersChange(
        time=T_6, growth_rate=-0.0022129, population_id=0),
    msprime.PopulationParametersChange(
        time=T_7, growth_rate=-0.0007438, population_id=0),
    msprime.PopulationParametersChange(
        time=T_8, growth_rate=-0.0016824, population_id=0),
    msprime.PopulationParametersChange(
        time=T_9, growth_rate=-0.0006301, population_id=0),
    msprime.PopulationParametersChange(
        time=T_10, growth_rate=-0.0005945, population_id=0),
    msprime.PopulationParametersChange(
        time=T_11, growth_rate=-0.0005306, population_id=0),
    msprime.PopulationParametersChange(
        time=T_12, growth_rate=-0.0000434, population_id=0),
    msprime.PopulationParametersChange(
        time=T_13, growth_rate=-0.0000, population_id=0),
    msprime.PopulationParametersChange(
        time=T_14, growth_rate=-0.0, population_id=0),

]
demographic_events1 = [

    #msprime.PopulationParametersChange(
     #   time=4, initial_size=, growth_rate=-0.095894024, population_id=0),
    msprime.PopulationParametersChange(
        time=4, initial_size=120, growth_rate=-0.24465639, population_id=0),
    msprime.PopulationParametersChange(
        time=7, growth_rate=-0.0560787, population_id=0),
    msprime.PopulationParametersChange(
        time=13, growth_rate=-0.1749704, population_id=0),
    msprime.PopulationParametersChange(
        time=19, growth_rate=-0.0675775, population_id=0),
    msprime.PopulationParametersChange(
        time=25, growth_rate=-0.0022129, population_id=0),
    msprime.PopulationParametersChange(
        time=155, growth_rate=-0.0007438, population_id=0),
    msprime.PopulationParametersChange(
        time=455, growth_rate=-0.0016824, population_id=0),
    msprime.PopulationParametersChange(
        time=655, growth_rate=-0.0006301, population_id=0),
    msprime.PopulationParametersChange(
        time=1755, growth_rate=-0.0005945, population_id=0),
    msprime.PopulationParametersChange(
        time=2355, growth_rate=-0.0005306, population_id=0),
    msprime.PopulationParametersChange(
        time=3355, growth_rate=-0.0000434, population_id=0),
    msprime.PopulationParametersChange(
        time=33155, growth_rate=-0.0000, population_id=0),
    msprime.PopulationParametersChange(
        time=933155, growth_rate=-0.0, population_id=0),

]

# Use the demography debugger to print out the demographic history
# that we have just described.
dd = msprime.DemographyDebugger(
    #  population_configurations=population_configurations,
    # migration_matrix=migration_matrix,
    demographic_events=demographic_events1)
dd.print_history()
recipe1 = pyslim.load("/home/iechevarriaz/SLiM/build/WholeChromosome90max3gen.trees").simplify()

print(recipe1.get_sequence_length())
print(recipe1.get_sample_size())
print(recipe1.get_num_nodes())
print(recipe1.get_time(471))

tree_sequence = msprime.simulate(   #sample_size= 12,
                                    Ne=0.25,
                                    length=1e8,
                                    from_ts=recipe1,
                                    recombination_rate=1e-8,
                                    #mutation_rate=2.5e-8,
                                    #population_configurations=population_configurations,
                                    #migration_matrix=migration_matrix,
                                    demographic_events=demographic_events1,
                                    random_seed= 248
                                )
mu_tree_seq = msprime.mutate(tree_sequence, rate= 1e-8, random_seed=786, keep=False)
#tree = mu_tree_seq.first()
#print(tree.draw(format="unicode"))

#for variant in mu_tree_seq.variants():
 #    print(variant.site.id, variant.site.position, variant.alleles, variant.genotypes, sep="\t")
print(mu_tree_seq.get_num_nodes(), mu_tree_seq.get_num_sites())
mu_tree_seq.dump("./HolsteinCattleTestcon90indWholeChrom90max.trees")

"""for tree in mu_tree_seq.trees():
    print("interval = ", tree.interval)
    print(tree.draw(format="unicode"))
"""
end = time.process_time()
diff = end - startFloat
CPUtime = str(diff)
corrCPU = str(end)

end2 = time.perf_counter()
diff2 = end2 - start2Float
WallTime = str(diff2)
corrWall = str(end2)
print("Elapsed CPU time: " + CPUtime)
print("Wall-Clock Elapsed time: " + WallTime)
print("CorroboroCPU" + corrCPU +"  Start: " + start)
print("CorroboroWall" + corrWall + "  Start2:  " + start2)