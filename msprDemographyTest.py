import msprime
import math
import pyslim
import time

#Clock Start
startFloat = time.process_time()
start = str(startFloat)
start2Float = time.perf_counter()
start2 = str(start2Float)

# If Times are provided in years, we convert into generations by replacing the number of years ago below.
generation_time = 5
st = 0
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

#population_configurations = [
  #  msprime.PopulationConfiguration(
  #      sample_size=10, initial_size=90),
  #  msprime.PopulationConfiguration(
  #      sample_size=10, initial_size=100),

#]

# Demographic history for Cattle
demographic_events = [

    msprime.PopulationParametersChange(
        # Here 'time' should be in generation notation, so if it is given in years, it should be replaced by T_0, T_1, ... , T_n; as shown above)
        # Growth rate is "per generation exponential growth rate": -alpha= [ln(initial_pop_size/next_stage_pop_size)/generation_span_in_years]
        # For example: ln(90/120)/3= -0.095894024
        time=1, initial_size=90, growth_rate=-0.095894024, population_id=0),
    msprime.PopulationParametersChange(
        time=4, growth_rate=-0.24465639, population_id=0),
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
    demographic_events=demographic_events)
dd.print_history()

tree_sequence = msprime.simulate(   sample_size= 2,
                                    #Ne=50,
                                    length=1e8,
                                    #from_ts=recipe1,
                                    start_time=T_0,
                                    #end_time=T_14,
                                    #__tmp_max_time=50,
                                    recombination_rate=1e-8,
                                    #mutation_rate=2.5e-8,
                                    #population_configurations=population_configurations,
                                    #migration_matrix=migration_matrix,
                                    demographic_events=demographic_events,
                                    random_seed=2478
                                )
#tree = tree_sequence.first()
mu_tree_seq = msprime.mutate(tree_sequence, rate= 1e-8, random_seed=45, keep=False)
#tree = mu_tree_seq.first()
#print(tree.draw(format="unicode"))
#for variant in tree_sequence.variants():
for variant in mu_tree_seq.variants():
     print(variant.site.id, variant.site.position, variant.alleles, variant.genotypes, sep="\t")

mu_tree_seq.dump("./DemographyTest.trees")

#for tree in tree_sequence.trees():
 #   print("interval = ", tree.interval)
  #  print(tree.draw(format="unicode"))

#Clock End
end = time.process_time()
diff = end - startFloat
CPUtime = str(diff)
checkCPU = str(end)

end2 = time.perf_counter()
diff2 = end2 - start2Float
WallTime = str(diff2)
checkWall = str(end2)
print("Elapsed CPU time: " + CPUtime)
print("Wall-Clock Elapsed time: " + WallTime)
print("CheckCPU" + checkCPU +"  Start: " + start)
print("CheckWall" + checkWall + "  Start2:  " + start2)

print(mu_tree_seq.get_num_sites(), mu_tree_seq.get_sample_size())
