import numpy as np


def multiplyArrayBy2(*args):
    """
    :brief Multiply the array by 2
    :param x: a np.ndarray
    """
    x = 2*np.array(args[0])
    print(x)
    return x
