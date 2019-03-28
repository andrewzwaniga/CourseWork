import sympy as sp

class Scalar(object):
    """Class defining a Scalar object. 

    A scalar is formally a (0,0) tensor.  

    Parameters: 
    -----------
    Pass parameters to the constructor that will be used to create attributes of a class 
    instance. Some parameters may or may not be attributes themselves. 

    index_dict : dict, string 
        A dictionary of strings representing the names of the variables to be used.
        For example, to create the Schwarzchild metric one should pass 
        {0:'t', 1:'r', 2:'theta', 3:'phi'}. 
    
    Attributes: 
    -----------
    Properties of a class instance that are calculated from constructor data. 

    self.elements : dict of float 
        Made to resemble the other classes used in this analysis,
        so that similar functions can be used on this even though 
        it has zero indices. We treat it as a dict where the key 
        is the empty string. 
    """

    def __init__(self, index_dict): 
        self.index_dict = index_dict 
        self.index_dim = len(self.index_dict) 
        self.elements = {'': 0.0}


