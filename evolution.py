'''
Handles how to perform all of the actual evolution.
'''
import random
import sys
from copy import copy
from util import diff_count


class Individual(object):
    '''
    An individual object used to combine gene fitness with genomes, as
    well methods for manipulating those genomes.
    '''
    def __init__(self, graph_length, input_length, output_length,
                  max_arity, function_list, **_):
        '''
        Create a new individual instance.

        Parameters:

        - ``graph_length``: The number of nodes in the CGP encoded graph.
        - ``input_length``: The number of input variables.
        - ``output_length``: The number of output variables.
        - ``max_arity``: The maximum arity used by any function.
        - ``function_list``: The list of functions a node can use.
        '''
        self.node_step = max_arity + 1
        self.input_length = input_length
        self.graph_length = graph_length
        self.function_list = function_list
        self.output_length = output_length
        self.genes = [self.random_gene(index) for index in
                      range(graph_length * self.node_step + output_length)]
        self.determine_active_nodes()
        # If memory problems arise, make this globally shared
        self.scratch = [None] * (graph_length + self.input_length)
        self.fitness = -sys.maxint

    def random_gene(self, index, invalid=None):
        '''
        Determines a random gene value given a gene index.  If optional
        ``invalid`` option is used, the returned value will only be ``invalid``
        if the gene has no other valid values.

        Parameters:

        - ``index``: The gene index who's value is being set.
        - ``invalid``: Value to avoid returning if possible
        '''
        node_number = index // self.node_step
        gene_number = index % self.node_step
        # If the gene is used to specify output locations
        if node_number >= self.graph_length:
            node_number = self.graph_length
            gene_number = -1
        # If the gene controls the function of a node
        if gene_number == 0:
            if len(self.function_list) == 1:
                return self.function_list[0]
            while True:
                choice = random.choice(self.function_list)
                if choice != invalid:
                    return choice
        # If the gene controls a connection / output location
        else:
            if node_number + self.input_length == 1:
                return -1
            while True:
                choice = random.randrange(-self.input_length, node_number)
                if choice != invalid:
                    return choice

    def copy(self):
        '''
        Return a copy of the individual.  Note that individuals are shallow
        copied except for their list of genes.
        '''
        # WARNING individuals are shallow copied except for things added here
        new = copy(self)
        new.genes = list(self.genes)
        return new

    def connections(self, node_index):
        '''
        Return the list of connections that a specified node has.

        Parameters

        - ``node_index``: The index of the node information is being requested
          for.  Note this is different from gene index.
        '''
        node_start = self.node_step * node_index
        return self.genes[node_start + 1: node_start + self.node_step]

    def determine_active_nodes(self):
        '''
        Determines which nodes are currently active and sets self.active
        to the sorted list of active genes.  Automatically called by gene
        manipulating member functions.
        '''
        self.active = set(self.genes[-self.output_length:])

        for node_index in reversed(range(self.graph_length)):
            if node_index in self.active:
                # add all of the connection genes for this node
                self.active.update(self.connections(node_index))
        self.active = sorted([acting for acting in self.active if acting >= 0])

    def evaluate(self, inputs):
        '''
        Given a list of inputs, return a list of outputs from executing
        this individual.

        Parameters:

        - ``inputs``: The list of input values for the individual to process.
        '''
        # Start by loading the input values into scratch
        # NOTE: Input locations are given as negative values
        self.scratch[-len(inputs):] = inputs[::-1]
        # Loop through the active genes in order
        for node_index in self.active:
            function = self.genes[node_index * self.node_step]
            args = [self.scratch[con] for con in self.connections(node_index)]
            # Apply the function to the inputs from scratch, saving results
            # back to the scratch
            self.scratch[node_index] = function(*args)
        # Extract outputs from the scratch space
        return [self.scratch[output]
                for output in self.genes[-self.output_length:]]

    def mutate(self, mutation_rate):
        '''
        Return a mutated version of this individual using the specified
        mutation rate.

        Parameters:

        - ``mutation_rate``: The probability that a specific gene will mutate.
        '''
        mutant = self.copy()
        for index in range(len(mutant.genes)):
            if random.random() < mutation_rate:
                mutant.genes[index] = mutant.random_gene(index,
                                                         mutant.genes[index])
        # Have the mutant recalculate its active genes
        mutant.determine_active_nodes()
        return mutant

    def one_active_mutation(self, _):
        '''
        Return a mutated version of this individual using the ``Single``
        mutation method.
        '''
        mutant = self.copy()
        while True:
            # Choose an index at random
            index = random.randrange(len(mutant.genes))
            # Get a new value for that gene
            newval = mutant.random_gene(index)
            # If that value is different than the current value
            if newval != mutant.genes[index]:
                mutant.genes[index] = newval
                # Determine if that gene was part of an active node
                node_number = index // self.node_step
                if (node_number >= self.graph_length or
                    node_number in self.active):
                    break
        # Have the mutant recalculate its active genes
        mutant.determine_active_nodes()
        return mutant

    def asym_phenotypic_difference(self, other):
        '''
        Determine how many genes would have to be mutated in order to make
        the ``other`` individual phenotypically identical to ``self``.

        Parameters:

        - ``other``: The individual to compare with.
        '''
        # Count the differences in the output locations
        count = diff_count(self.genes[-self.output_length:],
                           other.genes[-self.output_length:])
        # For each active node
        for node_index in self.active:
            index = node_index * self.node_step
            # Count the number of different connection genes
            count += diff_count(self.connections(node_index),
                                other.connections(node_index))
            # Include differences in the function gene
            count += (self.genes[index] !=
                      other.genes[index])
        return count

    def show_active(self):
        '''
        Prints the active portions of the individual in a somewhat readable
        way.
        '''
        for node_index in self.active:
            node_start = self.node_step * node_index
            print node_index, self.genes[node_start],
            print self.connections(node_index)
        print self.genes[-self.output_length:]

    def __lt__(self, other):
        '''
        Returns the result of self.fitness < other.fitness.
        '''
        return self.fitness < other.fitness

    def __le__(self, other):
        '''
        Returns the result of self.fitness <= other.fitness.
        '''
        return self.fitness <= other.fitness


def generate(config, output):
    '''
    An ``Individual`` generator that will yield a never ending supply of
    ``Individual`` objects that need to have their fitness set before the
    next ``Individual`` can be yielded.

    Parameters:

    - ``config``: A dictionary containing all configuration information
      required to generate initial individuals.  Should include values
      for:

      - All configuration information required to initialize an Individual.
      - ``mutation_rate``: The probably to use for mutation.
      - ``off_size``: The number of offspring to produce per generation.
      - ``output_length``: The number of output variables.
      - ``max_arity``: The maximum arity used by any function.
      - ``speed``: String specifying the way to handle duplicate
        individual creation, either ``normal'', ``skip'', ``accumulate``, or
        ``single``.
    - ``output``: Dictionary used to return information about evolution, will
      send out:

      - ``skipped``: The number of individuals skipped by ``Skip``.
      - ``estimated``: The estimated number of individuals that are skippable.
    '''
    output['skipped'] = 0
    output['estimated'] = 0
    if config['speed'] == 'single':
        Individual.mutate = Individual.one_active_mutation
    parent = Individual(**config)
    yield parent
    while True:
        mutants = [parent.mutate(config['mutation_rate'])
                   for _ in range(config['off_size'])]
        active = config['output_length'] + (len(parent.active) *
                                            (config['max_arity'] + 1))
        for index, mutant in enumerate(mutants):
            output['estimated'] += (1 - config['mutation_rate']) ** active
            prev = mutant
            if config['speed'] not in ['normal', 'single']:
                change = parent.asym_phenotypic_difference(mutant)
                if change == 0:
                    output['skipped'] += 1
                    if config['speed'] == 'skip':
                        continue
                    if config['speed'] == 'accumulate':
                        while change == 0:
                            prev = mutant
                            mutant = prev.mutate(config['mutation_rate'])
                            change = parent.asym_phenotypic_difference(mutant)
            yield mutant
            if config['speed'] == 'accumulate':
                # If the mutant is strickly worse, use the last equivalent
                mutants[index] = prev if mutant < parent else mutant
        best_child = max(mutants)
        if parent <= best_child:
            parent = best_child
