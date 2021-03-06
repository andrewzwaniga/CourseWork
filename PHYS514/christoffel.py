import sympy as sp

class Christoffel(object): 
    """Class defining a Christoffel object. 

    The Christoffel symbols are calculated from a metric as 

    \[
    '\Gamma^{\rho}_{\mu\nu} = \frac{1}{2}g^{\rho\lambda}(\partial_\mu g_{\nu\lambda} + \partial_\nu g_{\lambda\mu} - \partial_\lambda g_{\mu\nu})'
    \]

    Parameters: 
    -----------
    Pass parameters to the constructor that will be used to create attributes of a class 
    instance. Some parameters may or may not be attributes themselves.

    index_dict : dict, string 
        A dictionary of strings representing the names of the variables to be used. 
        For example, if working in spherical polar coordinates use 
        {0:'t', 1:'r', 2:'theta', 3:'phi'} 

    Attributes: 
    -----------
    Properties of a class instance that are calculated from constructor data. 

    self.elements : dict of dict of dict of float 
    """

    def __init__(self, index_dict): 
        self.index_dict = index_dict 
        self.index_dim = len(index_dict) 
        self.elements = {}
        for i in range(self.index_dim):
            self.elements[index_dict[i]] = {}
            for j in range(self.index_dim): 
                self.elements[index_dict[i]][index_dict[j]] = {}
                for k in range(self.index_dim):
                    self.elements[self.index_dict[i]][self.index_dict[j]][self.index_dict[k]] = 0.0

    def convert_to_shorthand(self): 
        """ A method that converts ``self.elements`` structure into a single
        dictionary, whereby $\Gamma^{t}_{tt}$ can be accessed with key ``['ttt']``. 
        Returns the conversion, but does not modify the class instance itself. 
        """

        shortcut = {}
        for i in range(self.index_dim):
            for j in range(self.index_dim):
                for k in range(self.index_dim):
                    s = '{a}{b}{c}'
                    s = s.format(a=self.index_dict[i],
                                 b=self.index_dict[j],
                                 c=self.index_dict[k])
                    shortcut[s] = self.elements[self.index_dict[i]][self.index_dict[j]][self.index_dict[k]]
        self.elements = {} # start fresh
        for key in shortcut:
            self.elements[key] = shortcut[key] 

