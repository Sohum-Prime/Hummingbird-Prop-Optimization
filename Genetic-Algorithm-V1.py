import random

# Define the fitness function to evaluate the quality of a solution


def fitness_function(solution):
    score = 0
    for i in range(len(solution)):
        if solution[i] == 1:
            score += 1
    return score

# Initialize the population with random solutions


def initialize_population(pop_size, solution_size):
    population = []
    for i in range(pop_size):
        individual = []
        for j in range(solution_size):
            individual.append(random.randint(0, 1))
        population.append(individual)
    return population

# Select the best individuals for reproduction


def selection(population, fitness_values, elite_size):
    elite = []
    for i in range(elite_size):
        max_index = fitness_values.index(max(fitness_values))
        elite.append(population[max_index])
        population.pop(max_index)
        fitness_values.pop(max_index)
    return elite

# Reproduce the next generation of solutions


def reproduce(elite, pop_size):
    next_population = []
    for i in range(pop_size):
        parent1 = random.choice(elite)
        parent2 = random.choice(elite)
        child = []
        for j in range(len(parent1)):
            random_num = random.uniform(0, 1)
            if random_num < 0.5:
                child.append(parent1[j])
            else:
                child.append(parent2[j])
        next_population.append(child)
    return next_population

# Introduce random mutations to the solutions


def mutate(population, mutation_rate):
    for i in range(len(population)):
        random_num = random.uniform(0, 1)
        if random_num < mutation_rate:
            mutation_index = random.randint(0, len(population[i])-1)
            population[i][mutation_index] = 1 - population[i][mutation_index]
    return population

# Main function to run the genetic algorithm


def genetic_algorithm(pop_size, solution_size, elite_size, mutation_rate, max_generations):
    population = initialize_population(pop_size, solution_size)
    for i in range(max_generations):
        fitness_values = [fitness_function(solution)
                          for solution in population]
        elite = selection(population, fitness_values, elite_size)
        population = reproduce(elite, pop_size)
        population = mutate(population, mutation_rate)
    best_solution = population[fitness_values.index(max(fitness_values))]
    return best_solution


# Test the genetic algorithm
result = genetic_algorithm(50, 100, 10, 0.01, 100)
print("Best solution:", result)
print("Fitness value:", fitness_function(result))
