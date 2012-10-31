'''
Defines each of the benchmark problems used as well as the function sets
for those problems.
'''
import operator
import itertools
import math


def nand(x, y):
    '''
    Simple Nand function for inclusion in function sets.
    '''
    return not (x and y)


def nor(x, y):
    '''
    Simple Nor function for inclusion in function sets.
    '''
    return not (x or y)


# Standard lists of operators for different problems to use
binary_operators = [operator.or_, operator.and_, nand, nor]


class Problem(object):
    '''
    Object used to store training input values and expected output values
    which are used in evaluating individuals.
    '''
    def __init__(self, problem_function, config):
        '''
        Create a new problem.

        Parameters:

        - ``problem_function``: The ground truth function that maps input
          values to output values.  Should be decorated with the
          ``problem_attributes`` decorator.
        - ``config``: A dictionary containing the configuration information
          required to fully initialize the problem.  Should include values
          for:

          - Any configuration information required to construct the problem
            range.
          - ``epsilon``: The amount of allowed error on each test.
        '''
        self.training = [(inputs, problem_function(inputs))
                         for inputs in problem_function.range(config)]
        self.epsilon = config['epsilon']

    def get_fitness(self, individual):
        '''
        Return the fitness of an individual as applied to this problem.

        Parameters:

        - ``individual``: The individual to be evaluated.
        '''
        score = 0
        for inputs, outputs in self.training:
            answers = individual.evaluate(inputs)
            # Finds the average number of outputs more than epsilon away from
            # the correct output
            score += (sum(float(abs(answer - output) > self.epsilon)
                          for answer, output in zip(answers, outputs))
                      / len(outputs))
        # Returns the percentage of correct answers
        return 1 - (score / float(len(self.training)))


def problem_attributes(range_method, operators, max_arity):
    '''
    Decorator that adds attributes to problem functions.

    Parameters

    - ``range_method``: The function used to generate training input values.
    - ``operators``: The list of valid operators on this problem.
    - ``max_arity``: The maximum arity of all operators.
    '''
    def wrapper(wraps):
        wraps.range = range_method
        wraps.operators = operators
        wraps.max_arity = max_arity
        return wraps
    return wrapper


def binary_range(config):
    '''
    Given a dictionary specifying the ``input_length``, returns all binary
    values of that length.
    '''
    return itertools.product((0, 1), repeat=config['input_length'])


def single_bit_set(config):
    '''
    Creates the list of all possible binary strings of specified length
    with exactly one set bit.  ``config`` should specify the ``input_length``.
    '''
    return [map(int, '1'.rjust(i + 1, '0').ljust(config['input_length'], '0'))
            for i in range(config['input_length'])]


@problem_attributes(binary_range, binary_operators, 2)
def even_parity(inputs):
    '''
    Return the even parity of a list of boolean values.
    '''
    return [(sum(inputs) + 1) % 2]


@problem_attributes(binary_range, binary_operators, 2)
def binary_multiply(inputs):
    '''
    Return the result of performing a binary multiplication of the first half
    of the inputs with the second half.  Will always have the same number of
    output bits as input bits.
    '''
    # convert the two binary numbers to integers
    joined = ''.join(map(str, inputs))
    middle = len(joined) / 2
    a, b = joined[:middle], joined[middle:]
    # multiply the two numbers and convert back to binary
    multiplied = bin(int(a, 2) * int(b, 2))[2:]
    # pad the result to have enough bits
    extended = multiplied.rjust(len(inputs), '0')
    return map(int, extended)


@problem_attributes(single_bit_set, binary_operators, 2)
def binary_encode(inputs):
    '''
    Returns the binary encoding of which input line contains a one.
    '''
    oneat = inputs.index(1)
    binary = bin(oneat)[2:]
    width = math.log(len(inputs), 2)
    return map(int, binary.zfill(int(width)))


@problem_attributes(binary_range, binary_operators, 2)
def binary_decode(inputs):
    '''
    Returns a 1 on the output line specified by the binary input index
    '''
    combined = ''.join(map(str, inputs))
    width = 2 ** len(inputs)
    base = [0] * width
    base[int(combined, 2)] = 1
    return base
