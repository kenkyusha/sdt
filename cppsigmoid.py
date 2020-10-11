import math

def cppsigmoid(x):

    sigma = 1./( 1 + math.exp(-x))
    return sigma

