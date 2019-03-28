import sympy as sp

class Ricci(object):
    """Class defining a Ricci object. 

    The Ricci tensor is a (0,2) tensor $R_{\mu\nu}$ given by the following trace of 
    the Riemann tensor: 
    
    \[
    R_{\mu\nu} = R^{\lambda}_{~\mu\lambda\nu}
    \]

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

    self.elements : dict of dict of float 
        The entries of the Ricci tensor, made to be accessible using a method that 
        resembles symbolic notation. That is, to get the $R_{tt}$ entry, 
        call ``self.elements['t']['t']`` 
    """

    def __init__(self, index_dict): 
        self.index_dict = index_dict 
        self.index_dim = len(self.index_dict) 
        self.elements = {}
        for i in range(self.index_dim):
            self.elements[index_dict[i]] = {}
            for j in range(self.index_dim):
                self.elements[index_dict[i]][index_dict[j]] = 0.0
    
    def convert_to_shorthand(self): 
        """ A method that converts ``self.elements`` structure into a single
        dictionary, whereby $R_{tt}$ can be accessed with key ``['tt']``. 
        Returns the conversion, but does not modify the class instance itself. 
        """

        shortcut = {}
        for i in range(self.index_dim):
            for j in range(self.index_dim):
                s = '{a}{b}'
                s = s.format(a=self.index_dict[i],
                             b=self.index_dict[j])
                shortcut[s] = self.elements[self.index_dict[i]][self.index_dict[j]]
        self.elements = {} # start fresh
        for key in shortcut:
            self.elements[key] = shortcut[key] 

