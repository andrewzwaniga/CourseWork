import sympy as sp

class Metric(object):
    """Class defining a Metric object. 

    A metric is a (0, 2) tensor whose components define the spacetime interval: 
    when the metric is denoted $g_{\mu\nu}$ then

    \[
    ds = \sqrt{g_{\mu\nu}\dot{x}^\mu\dot{x}^\nu}d\lambda
    \]
    
    where "dot" denotes differentiation with respect to the parameter $\lambda$ that
    labels points on the worldline. $ds^2$ is the spacetime interval between infinitesimally
    separated events. 

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
        The entries of the metric, made to be accessible using a method that 
        resembles symbolic notation. That is, to get the $g_{tt}$ entry, 
        call ``self.elements['t']['t']``. 
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
        dictionary, whereby $g_{tt}$ can be accessed with key ``['tt']``. 
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

