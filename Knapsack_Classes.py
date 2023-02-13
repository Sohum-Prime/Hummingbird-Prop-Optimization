"""
Class for knapsack content combos
"""

import random

class Thing:
    def __init__(name, value, weight):
        self.name = name
        self.value = value
        self.weight = weight

class Bag:
    possibilities = #list of Things
    
    def __init__(genome): #genome is a binary list
        self.num_items = len(genome)
        self.contents = [] #list of Things

        for i in range(num_items - 1):
            if genome[i] == 1:
                self.contents.append(possibilities[i])

    def eval_fitness(self, weight_limit):
        """
        Evaluates fitness of one knapsack object
        """

        weight_tot = 0
        value_tot  = 0
        for item in self.contents:
            weight_tot += item.weight
            value_tot  += item.value
            if weight_tot > weight_limit:
                return 0
        self.fitness = value_tot
        return value_tot

    def mutate(self, probability = 0.1):
        """
        Mutates the Bag
        """
        
        for i in range(self.num_items - 1):
            if random.random() < probability:
                self.genome[i] = abs(self.genome[i] - 1)

    def crossover(bag1, bag2):
        """
        Crosses over 2 bags (neither need to be self)

        args: 2 bags
        returns: 2 bags
        """

        cut_point = random.randint(1, length - 1)
        kid_genome1 = bag1.genome[0:cut_point] + bag2.genome[cut_point:]
        kid_genome2 = bag2.genome[0:cut_point] + bag1.genome[cut_point:]
        kid1 = Bag(kid_genome1)
        kid2 = Bag(kid_genome2)
        return kid1, kid2
        
