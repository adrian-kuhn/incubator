#!/usr/bin/env python
# -*- coding: utf-8 -*-

import random
import copy
import datetime
import os
import numpy
import arcpy
import logging

from scipy import spatial
from multiprocessing import Pool

log = logging.getLogger()


class Gene:
    """
    A gene represents a single parameter of an algorithm.
    """

    def __init__(self, name, pool, value=None):
        """
        Initialize a gene with name, gene pool and optional starting value.
        If no value provided, the gene will be initialized with a random value.

        :param name: The unique name of the gene
        :type name: String
        :param pool: The gene pool containing an interval of allowed values.
        :type pool: List of any values
        :param value: The value of the gene
        :type value: Any
        """
        self.name = name
        self.pool = pool
        self.value = value if value is not None else self.init()

    def __repr__(self):
        """
        String representation of this gene.

        :return: Gene as a String
        :rtype: String
        """
        return "{}: {}".format(self.name, round(self.value, 1) if isinstance(self.value, float) else self.value)

    def __eq__(self, other):
        """
        Genes are equal by name. Every gene in a genome needs a unique name to be distinguishable.

        :param other: The other gene to check with.
        :type other: :class:`Gene`
        :return: True, if this gene is identical with the other one. False otherwise
        :rtype: Boolean
        """
        return isinstance(other, self.__class__) and self.name == other.name

    def init(self):
        """
        Init a random value by selecting an index configured in the gene pool.

        :return: The newly initialized gene value.
        :rtype: Number
        """
        return self.pool[random.randint(0, len(self.pool) - 1)]

    def mutate(self):
        """
        Randomly mutate the value of this gene.
        """
        new_value = self.init()
        if new_value != self.value:
            self.value = new_value
        else:
            self.mutate()


class Chromosome:
    """
    A chromosome represents a collection of genes and therefore a collection of parameters.
    """

    def __init__(self, genes=None):
        """
        Initialize a chromosome with or without genes.

        :param genes: Optional list of genes.
        :type genes: List of :class:`Gene`
        """
        self.genes = genes

    @classmethod
    def random(cls, gene_pool):
        """
        Use this class function to init the chromosome with a gene pool dictionary.
        All the parameters will be initialized randomly.

        :param gene_pool: A dictionary with key value pairs containing genes names and possible values.
        :type gene_pool: Dictionary
        :return: The newly created chromosome with random gene values
        :rtype: :class:`Chromosome`
        """
        genes = []
        for name, pool in gene_pool.items():
            genes.append(Gene(name, pool))
        return cls(genes)

    def __repr__(self):
        """
        String representation of this chromosome

        :return: Chromosome as a String
        :rtype: String
        """
        return " | ".join(str(gene) for gene in self.genes)

    def __eq__(self, other):
        """
        Equality check for two chromosomes. Chromosomes are identical by the genes.

        :param other: The other chromosome to check with.
        :type other: :class:`incubator.Chromosome`
        :return: True, if both chromosomes are identical. False otherwise
        :rtype: Boolean
        """
        return isinstance(other, self.__class__) and self.genes == other.genes

    def mutate(self, number_of_genes=3):
        """
        Randomly mutate genes (one parameter) in the chromosome (in the set of parameters)

        :param number_of_genes: Define the number of genes being mutated randomly. Default is 3.
        :type number_of_genes: Integer
        """
        for i in range(0, number_of_genes):
            self.genes[random.randint(0, len(self.genes) - 1)].mutate()

    def crossover(self, other, erroneous=True):
        """
        Crossover function to mix this chromosome with another one by dividing the genome in two half.
        The child chromosome is built by concatenating half of the father's and half of the mother's genome.

        :param other: The other chromosome to mix this chromosome with
        :type other: :class:`Chromosome`
        :param erroneous: Define if this genome crossover function is erroneous and can have mutations. Default is True.
        :type erroneous: Boolean.
        :return: New chromosome build with parameters from this and the other chromosome
        :rtype: :class:`Chromosome`
        """
        if self == other:
            parent_1 = self.genes[:len(self.genes) // 2]
            parent_2 = other.genes[len(other.genes) // 2:]
            chromosome = Chromosome(copy.deepcopy(parent_1) + copy.deepcopy(parent_2))
            if erroneous:
                chromosome.mutate()
            return chromosome


class Individual:
    """
    Abstract base class for a genetic algorithm with a chromosome (list of parameters).
    """

    def __init__(self, chromosome=None, gene_pool=None, assignments=None):
        """
        Constructor can be used with a predefined chromosome.

        :param chromosome: An optional existing chromosome to init this individual
        :type chromosome: :class:`Chromosome`
        :param gene_pool: Static Dictionary with intervals for every gene.
        :type gene_pool: Dictionary
        :param assignments: List of assignments to solve
        :type assignments: List of :class:`Assignment`
        """
        if chromosome is None and gene_pool is None:
            raise AttributeError("Genetic algorithms need a either a chromosome or a gene pool for initializing")

        self.chromosome = chromosome if chromosome is not None else Chromosome.random(gene_pool)
        self.assignments = assignments
        self.fitness = 0
        self.true_positives = 0
        self.false_positives = 0
        self.false_negatives = 0
        self.distance = 0
        self.fitness_weighted = 0

    def __eq__(self, other):
        """
        Equality check for two individuals. Individuals are identical by the chromosome.

        :param other: The other individual to check
        :type other: :class:`incubator.Individual`
        :return: True, if both individuals are identical. False otherwise
        :rtype: Boolean
        """
        return isinstance(other, self.__class__) and self.chromosome == other.chromosome

    def __repr__(self):
        """
        String representation of this individual.

        :return: Algorithm as a String
        :rtype: String
        """
        return "{} with chromosome {} and avg fitness *{}*, weighted fitness !{}! (Distance={}, TP={}, FP={}, FN={})"\
            .format(self.__class__.__name__, self.chromosome, self.fitness, self.fitness_weighted, self.distance,
                    self.true_positives, self.false_positives, self.false_negatives)

    def __gt__(self, other):
        """
        Makes algorithms sortable. An algorithm is greater than another one, if it has a smaller fitness value.

        :param other: The other individual to compare this algorithm with
        :type other: :class:`GeneticAlgorithm`
        :return: True, if this algorithm is fitter than the other one, false otherwise.
        :rtype: Boolean
        """
        return self.fitness > other.fitness

    def read_chromosome(self):
        """
        Abstract method to read the specific chromosome values. Must be implemented in concrete algorithm,
        because every algorithm has another chromosome.

        :raises: NotImplemented Exception, because it is abstract
        """
        raise NotImplemented("Must be implemented in concrete algorithm")

    def solve_assignment(self, **raw_data):
        """
        Abstract method to solve arbitrary assignment. Must be implemented in concrete algorithm.

        :param raw_data: Key worded arguments according settings needed to solve the assignment.
        :type raw_data: Any
        :raises: NotImplemented Exception, because it is abstract
        """
        raise NotImplemented("Must be implemented in concrete algorithm")

    def mate(self, other):
        """
        Abstract method to mate this individual with another one. Must be implemented in concrete algorithm.

        :param other: The other individual
        :type other: :class:`incubator.Individual`
        :raises: NotImplemented Exception, because it is abstract
        """
        raise NotImplemented("Must be implemented in concrete algorithm")

    def run(self):
        """
        Solve all assignments stored in this individual and calculate F1-Score as fitness parameter.
        """
        self.read_chromosome()
        if self.fitness == 0:
            for assignment in self.assignments:
                assignment.result = self.solve_assignment(**assignment.raw_data)
                assignment.fit()
                self.fitness += assignment.result.fitness if assignment.result is not None else 0
                self.distance += assignment.result.distance if assignment.result is not None else 0
                self.true_positives += assignment.result.true_positives if assignment.result is not None else 0
                self.false_positives += assignment.result.false_positives if assignment.result is not None else 0
                self.false_negatives += assignment.result.false_negatives if assignment.result is not None else 0

            self.fitness = self.fitness / len(self.assignments)
            ua = self.true_positives / (self.true_positives + self.false_positives)
            pa = self.true_positives / (self.true_positives + self.false_negatives)
            if ua + pa > 0:
                self.fitness_weighted = 2 * ((ua * pa) / (ua + pa))


class Generation:
    """
    Collection of multiple individuals
    """
    def __init__(self, individuals=None):
        """
        Initialize the generation with an optional list of individuals.

        :param individuals: The concrete individuals (the algorithms solving the problem)
        :type individuals: List of class:`Individual`
        """
        self.individuals = individuals or []

    def __repr__(self):
        """
        String representation of this generation.

        :return: Generation as a String
        :rtype: String
        """
        return "\n".join(str(individual) for individual in self.individuals)

    def run_parallel(self, pool):
        """
        Solve all assignments within the individuals of this generation with multiple processes.

        :param pool: Processing pool from the python multiprocessing library.
        :type pool: Python multiprocessing pool
        """
        self.individuals = pool.map(self._run, self.individuals)

    def run_serial(self):
        """
        Solve all assignments within the individuals of this generation with a single process.
        """
        for individual in self.individuals:
            individual.run()

    @staticmethod
    def _run(individual):
        """
        Static function to solve the assignment in a parallel worker instance.

        :param individual: The individual
        :type individual: :class:`Individual`
        :return: The individual containing the results after solving the assignments.
        :rtype: :class:`Individual`
        """
        individual.run()
        return individual

    def evolve(self, number_of_parents, number_of_offsprings):
        """
        Select best performing individuals in this generation, mate them and build the next generation.

        :param number_of_parents: Defines the number of best performing parents being selected.
        :type number_of_parents: Integer
        :param number_of_offsprings: Defines the number of offsprings generated during the mating process.
        :type number_of_offsprings: Integer
        :return: The newly created generation containing the best performing parents and their offsprings.
        :rtype: :class:`Generation`
        """
        parents = self.selection(number_of_parents)
        offsprings = self.mating(parents, number_of_offsprings)
        return Generation(parents + offsprings)

    def selection(self, number_of_parents):
        """
        Select the best performing parents by sorting the individuals according the fitness parameter.

        :param number_of_parents: Defines the number of best performing parents being selected.
        :type number_of_parents: Integer
        :return: List of best performing individuals
        :type: List of class:`Individual`
        """
        return sorted(copy.deepcopy(self.individuals), reverse=True)[:number_of_parents]

    @staticmethod
    def mating(parents, number_of_offsprings):
        """
        A list of parents are mated by chromosome crossover and building offsprings.

        :param parents: The list of parents to mate
        :type parents: List of class:`Individual`
        :param number_of_offsprings: Defines the number of offsprings generated during the mating process.
        :type number_of_offsprings: Integer
        :return: List of newly created Individuals
        :rtype: List of class:`Individual`
        """
        number_of_parents = len(parents)
        offsprings = []
        for i in range(0, number_of_offsprings):
            offsprings.append(parents[i % number_of_parents].mate(parents[(i + 1) % number_of_parents]))
        return offsprings


class Incubator:
    """
    Main class containing a list of generations and controlling the evolutionary process.
    """
    def __init__(self, individual, assignments, number_of_generations, number_of_individuals, parallel=True):
        """
        Initialize the Incubator with settings

        :param individual: The concrete class used in the incubator (e.g. WatershedAdaptive)
        :type individual: class:`Individual`
        :param assignments: A list of assignments feed to the individuals in the incubator.
        :type assignments: class:`Assignment`
        :param number_of_generations: The maximum number of generations calculated in the incubator.
        :type number_of_individuals: Integer
        :param number_of_individuals: The number if individuals per generation used in the incubator.
        :type number_of_individuals: Integer
        :param parallel: Indication to solve the assignments serially or in parallel. Default is parallel (True)
        :type parallel: Boolean
        """
        self.assignments = assignments
        self.number_of_generations = number_of_generations
        self.number_of_individuals = number_of_individuals
        self.number_of_parents = number_of_individuals // 2
        self.number_of_offsprings = number_of_individuals // 2

        self.generations = []
        self.init_first_generation(individual)

        # Create a multiprocessing Pool
        self.parallel = parallel
        if self.parallel:
            self.pool = Pool(os.cpu_count())

    def init_first_generation(self, individual):
        """
        Initializes the incubator with a list of individuals in the first generation.
        The individuals will have random parameters.

        :param individual: Class instance containing a concrete representation of an individual (the algorithm).
        :type individual: class:`Individual`
        """
        individuals = []
        for i in range(0, self.number_of_individuals):
            individuals.append(individual(None, self.assignments))

        self.generations.append(Generation(individuals))

    def breed(self):
        """
        Looping over a specified number of generations an let the individuals solve the assignments.
        Select, mate and mutate the best performing individuals of the current generation and build the next generation.
        """
        for i in range(0, self.number_of_generations):
            start = datetime.datetime.now()
            log.info("Generation {}".format(i))
            if self.parallel:
                self.generations[i].run_parallel(self.pool)
            else:
                self.generations[i].run_serial()

            log.info(self.generations[i])
            next_generation = self.generations[i].evolve(self.number_of_parents, self.number_of_offsprings)
            self.generations.append(next_generation)
            end = datetime.datetime.now()
            log.info("Calculation time: {}".format(end - start))


class Assignment:
    """
    Represents an assignment containing the required raw data and target data set (the reference)
    """
    def __init__(self, target, **raw_data):
        """
        Initialize the assignment with a target object (the reference) and
        key worded arguments needed to solve the assignment.

        :param target: The target object (the reference data set)
        :type target: class:`TreeTops`
        :param raw_data: Arbitrary list of data needed so solve the assignment.
        :type raw_data: Any
        """
        self.target = target
        self.raw_data = raw_data
        self.result = None

    def fit(self):
        """
        Calculate the fitness parameter.
        """
        if self.result:
            self.result.similarity(self.target)


class TreeTops:
    """
    Representation of a point cloud with tree crown points (x, y, z).
    """
    def __init__(self, points):
        self.points = points
        """Numpy array with 3D point cloud"""

        self.x = points[:, 0]
        """First column represents values on x coordinate"""

        self.y = points[:, 1]
        """Second column represents values on y coordinate"""

        self.z = points[:, 2]
        """Third column represents values on z coordinate"""

        self.true_positives = 0
        self.false_positives = 0
        self.false_negatives = 0
        self.fitness = 0
        self.distance = 0

    @classmethod
    def from_feature_class(cls, path):
        """
        Load a point cloud from ESRI file geo database point feature class

        :param path: Path to the feature class on the file system
        :type path: String
        :return: A TreePoint object containing the tree crown points according the point feature class.
        :rtype: class:`TreeTops`
        """
        fields = ['SHAPE@XYZ']
        points = []
        with arcpy.da.SearchCursor(path, fields) as cursor:
            for row in cursor:
                points.append(row[0])

        return cls(numpy.array(points))

    def similarity(self, other, threshold=0.5):
        """
        Calculate similarity measurement of two tree point clouds.

        Measures the minimum euclidean distance for every point in target to the next point in test point cloud and vice versa.
        If the distance is equal to 0, two points match (TP = true positives).
        If the distance is greater than 0, the points doesn't match. We have a false positive (FP) or a false negative (FN)
        depending of the direction the euclidean distance was calculated.

        With a confusion matrix the accuracy index will be calculated according the equation TP / (TP + FP + FN)
        To get a weighted accuracy value, not the amount of FP or FN will be used, but the sum of the euclidean distances.

        :param other: The other point cloud to compare this point cloud with
        :type other: :class:`TreeTops`
        :param threshold: The threshold in meter the treat two points as equal. Default value is 0.0
        :type threshold: Float
        :return: Accuracy measurement in percent
        :rtype: Float
        """
        if len(self.points) != 0:
            min_dist_false_pos = numpy.min(numpy.around(spatial.distance.cdist(self.points, other.points), 2), axis=1)
            min_dist_false_neg = numpy.min(numpy.around(spatial.distance.cdist(other.points, self.points), 2), axis=1)

            self.true_positives = len(min_dist_false_pos[numpy.where(min_dist_false_pos <= threshold)])
            self.false_positives = len(min_dist_false_pos[numpy.where(min_dist_false_pos > threshold)])
            self.false_negatives = len(min_dist_false_neg[numpy.where(min_dist_false_neg > threshold)])

            # distance fitness
            self.distance = numpy.sum(min_dist_false_pos) + numpy.sum(min_dist_false_neg)

            # F-score
            ua = self.true_positives / (self.true_positives + self.false_positives)
            pa = self.true_positives / (self.true_positives + self.false_negatives)
            if ua + pa == 0:
                self.fitness = 0
            else:
                self.fitness = 2 * ((ua * pa) / (ua + pa))
        else:
            self.false_negatives = len(self.points)
            self.fitness = 0
