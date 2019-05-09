'''
Created on Jun 30, 2015

@author: javi

Runner Script for using the simuPOP library. Used to generate simulated recombinant genomes.
'''

import simuPOP as sp
import Simulation.SNPParser as snp
import ChrNameSorter as cns
import random as rand
import string
import os
import Simulation.Recombinator as rec
import functools as fc

def generateAncestors(n, num_chrs, snp_per_chr):
    '''
    generates some dummy ancestors. Super fake. 
    '''
    
    ancestors = {}
    nucs = ['A', 'T']
    for x in range(n):
        temp = ancestors[rec.formatName(str(x))] = {}
        for n in range(num_chrs):
            temp[chrName(n + 1)] = [rand.choice(nucs) for x in range(snp_per_chr)]
    
    pos_tree = {chr:[x for x in range(0,snp_per_chr * 50, 50)] for chr in temp}
    
    return ancestors, pos_tree

#Define digit mapping


def chrName(n):
    romanNumeralMap = (
                   ('X',  10),
                   ('IX', 9),
                   ('V',  5),
                   ('IV', 4),
                   ('I',  1))
    
    """convert integer to Roman numeral"""
    result = ""
    for numeral, integer in romanNumeralMap:
        while n >= integer:
            result += numeral
            n -= integer
    return "TGME49_CHR" + result

def parseSNPInput(file_path, num):
    data_complex = snp.parse(file_path)
    return snp.select(data_complex, num)

def createRateMatrix(dimension, rate):
    '''
    creates a dummy rate matrix with uniform rates for
    everything
    '''
    self_rate = 1 - rate * (dimension - 1)
    if self_rate < 0: raise Exception('createRateMatrix: rate is too high')
    
    base = [[rate] * dimension for x in range(dimension)]
    for x in range(dimension):
        base[x][x] = 0
    
    return base

def addNames(name):
    chars = sorted(set("".join([str(int(x)) for x in name])))
    return round(float("".join(chars)),2)


def generateID(length):
    '''
    generates an unique ID of set length
    '''
    return ''.join(rand.sample(string.ascii_uppercase, 3))
    
def addLists(x, y):
    return x + y    
    
def simulateBySimuPOP():
    #starting variables
    directory = '/data/new/javi/toxo/simulations4/'
    input_path = 'Toxo20.txt'
    output_path = 'SimulatedToxo.txt'
    
    input_path = directory + input_path
    output_path = directory + output_path
    parents_path = directory + '/parents.txt'
    pedigree_path = directory + 'pedigree.txt'

    
    
    number_of_ancestors = 3
    expansion_pop_size = 15
    offsprings_sampled = number_of_ancestors + expansion_pop_size
    gen = 3
    translate_mode = 'toxoplasma'
    structure_mode = 'simupop'
    
    #parsing input
    init_info = parseSNPInput(input_path, number_of_ancestors)
    ancestral_genomes = init_info[0]
    ancestor_names = ancestral_genomes.keys()
    loci_positions = init_info[1]
    chromosome_names = sorted(loci_positions.keys(), key = lambda x: cns.getValue(x, translate_mode))
    list_of_loci = [len(loci_positions[chr]) for chr in chromosome_names]
    lociPos = fc.reduce(lambda x, y: x + y, [loci_positions[x] for x in chromosome_names])
    
    sp.turnOnDebug(code = "DBG_GENERAL")
    
    #initializing
    print('Initializaing Population')
    population = sp.Population(size=[number_of_ancestors], loci=list_of_loci, ancGen = 5, lociPos = lociPos, \
                               chromNames = chromosome_names, lociNames = [], alleleNames = ['A','T','G','C'],\
                               infoFields=['name', 'ind_id', 'father_id', 'mother_id'])
    
    for individual, sample, ind_id in zip(population.individuals(), ancestral_genomes, range(len(ancestral_genomes))):
        individual.setInfo(ancestor_names.index(sample), 'name')
        individual.setInfo(ind_id, 'ind_id')
        for ind, chr in enumerate(chromosome_names):
            individual.setGenotype(ancestral_genomes[sample][chr], chroms=[ind])
    
    #Alternating rounds of recombination with clonal expansion. Clonal expansion gives + 2.
    #Mutation prior to each round
    
    simulator = sp.Simulator(population)
    rate_matrix = createRateMatrix(len(ancestor_names), 0.0002) #10,000 times the mutation rate.
    id_tagger = sp.IdTagger()
    ped_tagger = sp.PedigreeTagger(output='>>' + pedigree_path, outputFields=['name', 'ind_id'])
    inherit_tagger = sp.InheritTagger(infoFields = 'name')
    
    initOps1 = [sp.PyExec('print("Starting random selection")'), ped_tagger]
    initOps2 = [sp.PyExec('print("Starting random mating")'), ped_tagger]
    preOps1 = [sp.MatrixMutator(rate = rate_matrix)
               ]
    preOps2 = [sp.InitSex(sex=[sp.MALE, sp.FEMALE])]

    matingScheme1 = sp.RandomSelection(ops = [
                                              sp.CloneGenoTransmitter(),
                                              inherit_tagger,
                                              id_tagger,
                                              ped_tagger
                                              ],
                                              subPopSize = expansion_pop_size)
    matingScheme2 = sp.RandomMating(ops = [
                                           sp.Recombinator(intensity = 0.01 / 105000, convMode=(sp.GEOMETRIC_DISTRIBUTION, 0.001,0.01)), #10x normal
                                           sp.PyTagger(func = addNames),
                                           id_tagger,
                                           ped_tagger
                                           ],
                                 subPopSize = expansion_pop_size)
    
    postOps = []
    finalOps = []
    
    print('Starting Evolution Cycles.')
    
    try:
        os.remove(pedigree_path)
    except:pass
    
    simulator.evolve(initOps = [id_tagger, ped_tagger], matingScheme = sp.CloneMating(ops = [sp.CloneGenoTransmitter(), ped_tagger, id_tagger, inherit_tagger]), gen = 1)
    
    for x in range(gen):
        simulator.evolve(initOps = initOps1, preOps = preOps1, matingScheme = matingScheme1, postOps = postOps, finalOps = finalOps, gen = 1)
        simulator.evolve(initOps = initOps2, preOps = preOps2, matingScheme = matingScheme2, postOps = postOps, finalOps = finalOps, gen = 1)
        
    offsprings = {''.join([str(int(x.info('name'))),generateID(3), str(int(x.info('ind_id')))]): x.info('ind_id') for x in simulator.population(0).individuals()}
    sampled_ones = rand.sample(offsprings.keys(), offsprings_sampled)

    
    #reorganizes the offspring genome. Extract info by chr.
    offspring_genomes = {name: {} for name in sampled_ones}
    for name in sampled_ones:
        for ind, chr in enumerate(chromosome_names):
            offspring_genomes[name][chr] = simulator.population(0).indByID(offsprings[name], idField='ind_id').genotype(ploidy=0,chroms=[ind])

    offspring_genomes.update(ancestral_genomes)
    
    print('Parent Guide:')
    for ind, id in enumerate(ancestor_names):
        print(" : ".join([str(ind), str(id)]))
    
    print('Complete. Generating Output.')
    
    with open(parents_path, 'w') as parent_output:
        parent_output.write('Parent Guide:\n')
        for ind, id in enumerate(ancestor_names):
            parent_output.write(" : ".join([str(ind), str(id)]) + '\n')
    
    #output
    offspring_genomes = snp.restructure((offspring_genomes, loci_positions), structure_mode)
    snp.outputGriggFormat(offspring_genomes, output_path)
    print('Simulation Complete.')
    
    
def simulateByRecombinator():
    
    #starting variables
    directory = '/data/new/javi/toxo/simulations5/'
    input_path = 'Toxo20.txt'
    output_path = 'SimulatedToxo.txt'
    
    input_path = directory + input_path
    output_path = directory + output_path
    events_path = directory + 'events.txt'

    number_of_chrs = 14
    snps_per_chr = 30000
    expansion_pop_size = 15
    translate_mode = 'toxoplasma'
    structure_mode = 'generated'
    
    #parsing input
    init_info = generateAncestors(4, number_of_chrs, snps_per_chr)
    population = init_info[0]
    ancestor_names = population.keys()
    pos_tree = init_info[1]
    events = []
    
    #recombine
    print('Starting recombination cycles...')
    cycle_result = rec.cycle(population, pos_tree, expansion_pop_size)
    population = cycle_result[0]
    events.append(cycle_result[1])

    
    print('Complete. Generating Output.')
    
    with open(events_path, 'w') as events_output:
        for section in events:
            for item in section:
                events_output.write("\t".join(list(item)) + '\n')
    
    #output
    offspring_genomes = snp.restructure((population, pos_tree), structure_mode)
    snp.outputGriggFormat(offspring_genomes, output_path)
    print('Simulation Complete.')

        
if __name__ == '__main__':
    print('Simulating by Recombinator')
    simulateByRecombinator()
#     print('Simulatioin by simuPOP')
#     simulateBySimuPOP()
    
    
    
    