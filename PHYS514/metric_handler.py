""" Metric Handler 

Take as input an arbitrary metric with two lower indices and arbitrary dimension. 
From the metric, calculate the 
(1) Christoffel coefficients 
(2) Riemann tensor 
(3) Ricci tensor, Ricci scalar, and Einstein Tensor 

"""

import sympy as sp

christoffel_latex = '\Gamma^{\rho}_{\mu\nu} = \frac{1}{2}g^{\rho\lambda}(\partial_\mu g_{\nu\lambda} + \partial_\nu g_{\lambda\mu} - \partial_\lambda g_{\mu\nu})' 

riemannT_latex = 'R^{\rho}_{~\sigma\mu\nu} = \partial_\mu\Gamma^{\rho}_{\nu\sigma} - \partial_\nu\Gamma^{\rho}_{\mu\sigma} + \Gamma^{\rho}_{\lambda\mu}\Gamma^{\rho}_{\nu\sigma} - \Gamma^{\rho}_{\lambda\nu}\Gamma^{\rho}_{\mu\sigma}'

ricciT_latex = 'R_{\mu\nu} = R^{\lambda}_{~\mu\lambda\nu}'

ricciS_latex = 'R = g^{\mu\nu}R_{\mu\nu}'

einsteinT = 'G_{\mu\nu} = R_{\mu\nu} - \frac{1}{2}Rg_{\mu\nu}' 
                               
class Metric():
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
        return shortcut 

class Christoffel(): 
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
        return shortcut 

class Riemann(): 
    """Class defining a Riemann object. 

    The Riemann tensor is calculated from Christoffel symbols as 

    \[
    'R^{\rho}_{~\sigma\mu\nu} = \partial_\mu\Gamma^{\rho}_{\nu\sigma} - \partial_\nu\Gamma^{\rho}_{\mu\sigma} + \Gamma^{\rho}_{\lambda\mu}\Gamma^{\rho}_{\nu\sigma} - \Gamma^{\rho}_{\lambda\nu}\Gamma^{\rho}_{\mu\sigma}'
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

    self.elements : dict of dict of dict of dict of float
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
                    self.elements[index_dict[i]][self.index_dict[j]][self.index_dict[k]] = {}
                    for l in range(self.index_dim):
                        self.elements[self.index_dict[i]][self.index_dict[j]][self.index_dict[k]][self.index_dict[l]] = 0.0

    def convert_to_shorthand(self): 
        """ A method that converts ``self.elements`` structure into a single
        dictionary, whereby $\R^{t}_{ttt}$ can be accessed with key ``['tttt']``. 
        Returns the conversion, but does not modify the class instance itself. 
        """

        shortcut = {}
        for i in range(self.index_dim):
            for j in range(self.index_dim):
                for k in range(self.index_dim):
                    for l in range(self.index_dim):
                        s = '{a}{b}{c}{d}'
                        s = s.format(a=self.index_dict[i],
                                     b=self.index_dict[j],
                                     c=self.index_dict[k],
                                     d=self.index_dict[l])
                        shortcut[s] = self.elements[self.index_dict[i]][self.index_dict[j]][self.index_dict[k]][self.index_dict[l]]
        return shortcut 

    
class MetricAnalysis():
    """Class with functions that will operate on a given metric
    to produce the Christoffel symbols, the Riemann tensor, the Ricci
    tensor, Ricci scalar, and Einstein tensor. 
    """
    
    def __init__(self):
        self.SCHWARZCHILD = {0:'t', 1:'r', 2:'theta', 3:'phi'}
        self.S2 = {0:'theta', 1:'phi'}
        self.DESITTER = {0:'t', 1:'x', 2:'y', 3:'z'} 
        self.KERR = {0:'t', 1:'r', 2:'theta', 3:'phi'}

    def init_schwarzchild(self):
        """Initialize Schwarzchild metric.
        """
        t = sp.symbols(self.SCHWARZCHILD[0]) 
        r = sp.symbols(self.SCHWARZCHILD[1])
        theta = sp.symbols(self.SCHWARZCHILD[2])
        phi = sp.symbols(self.SCHWARZCHILD[3])
        R = sp.symbols('R')

        g = Metric(index_dict=self.SCHWARZCHILD)
        g = g.convert_to_shorthand()
        g['tt'] = 1 - R/r
        g['rr'] = (1 - R/r)**(-1)
        g['thetatheta'] = r**2
        g['phiphi'] = r**2*(sp.sin(theta))**2

        return g 

    def init_S2(self):
        """Initialize the metric on the 2-sphere, S2. 
        """
        theta = sp.symbols(self.S2[0])
        phi = sp.symbols(self.S2[1])
        g = Metric(index_dict=self.S2)
        g = g.convert_to_shorthand()
        g['thetatheta'] = 1 
        g['phiphi'] = (sp.sin(theta))**2

        return g 

    def init_deSitter(self):
        """Initialize the metric on 4D de Sitter spacetime.
        """
        t = sp.symbols(self.DESITTER[0])
        x = sp.symbols(self.DESITTER[1])
        y = sp.symbols(self.DESITTER[2])
        z = sp.symbols(self.DESITTER[3])
        H = sp.symbols('H') 

        g = Metric(index_dict=self.DESITTER)
        g = g.convert_to_shorthand()
        g['tt'] = -1 
        g['xx'] = sp.exp(2*H*t)
        g['yy'] = g['xx'] 
        g['zz'] = g['xx'] 

        return g 

    def init_kerr(self):
        """Initialize the Kerr metric on 4D spacetime.
        """
        
        t = sp.symbols(self.KERR[0])
        r = sp.symbols(self.KERR[1])
        theta = sp.symbols(self.KERR[2])
        phi = sp.symbols(self.KERR[3])
        G = sp.symbols('G')
        M = sp.symbols('M')
        rho = sp.symbols('rho')
        a = sp.symbols('a')
        Delta = sp.symbols('Delta')
        
        """
        Notes:
        ------
        a = J/M where J is the Komar angular momentum (see 6.48 
        in Carroll). 
        """

        Delta_expr = r**2 - 2*G*M*r + a**2 
        rho2_expr = r**2 + a**2*(sp.cos(theta))**2

        g = Metric(index_dict=self.KERR)
        g = g.convert_to_shorthand()

        g['tt'] = -(1-2*G*M*r/rho**2)
        g['tt'] = g['tt'].subs(rho**2, rho2_expr)
        
        g['tphi'] = -2*G*M*a*r*(sp.sin(theta))**2/rho**2
        g['tphi'] = g['tphi'].subs(rho**2, rho2_expr)
        g['phit'] = g['tphi'] 
        
        g['rr'] = rho**2/Delta
        g['rr'] = g['rr'].subs(rho**2, rho2_expr)
        g['rr'] = g['rr'].subs(Delta, Delta_expr)
        
        g['thetatheta'] = rho**2
        g['thetatheta'] = g['thetatheta'].subs(rho**2, rho2_expr)
                
        g['phiphi'] = ((sp.sin(theta))**2/rho**2)*((r**2 + a**2)**2 - a**2*Delta*(sp.sin(theta))**2)
        g['phiphi'] = g['phiphi'].subs(rho**2, rho2_expr)
        g['phiphi'] = g['phiphi'].subs(Delta, Delta_expr)
        
        return g 

    def invert_metric(self, g, index_dict): 
        """ Invert a given metric g and return the inverted metric 
        as a shorthand dictionary of sympy.core.symbols.Symbol objects.
        """

        import numpy as np 

        n = int(np.sqrt(len(g))) # casting will not chop, we expect whole number anyways
        #g_matrix = sp.zeros(n,n)
        g_matrix = np.zeros([n,n])
        for i in range(n): 
            for j in range(n):
                m = index_dict[i]
                n = index_dict[j] 
                mn = m+n
                g_matrix[i][j] = g[mn] 
        
        # for my next trick, it seems unavoidable to NOT to this manually...
        g_matrix = sp.Matrix(g_matrix[0],
                             g_matrix[1],
                             g_matrix[2],
                             g_matrix[3])
        #g_matrix_inv = np.linalg.inv(g_matrix)
        g_inv = g_matrix**(-1)        
        g_inv = {}
        for i in range(n):
            for g in range(n): 
                m = index_dict[i]
                n = index_dict[j] 
                mn = m+n
                g_inv[mn] = g_matrix_inv[i][j]  

    def calculate_Christoffel(self, g, index_dict, simplify): 
        """Calculate the Christoffel symbols for a given metric ``g``. 
        
        Parmeters:
        ---------- 
        g : dict, sympy.core.symbol.Symbol 
            Construct a blank metric first with the Metric class. 
            Then use sympy to create the desired metric. 
            Lastly, pass the shorthand dictionary to this function. 

        index_dict : dict, string 
            The names of the indices as strings. 

        simplify : bool 
            If True, simplify with sp.simplify(). Else do nothing. 

        Returns: 
        -------- 
        Gamma : dict, sympy.core.symbol.Symbol 
            Through this calculation, ``Gamma`` is implicitly converted to a sympy.symbols() object. 
            At first it is initialized as a blank Christoffel() object then converted to shorthand. 
        """

        Gamma = Christoffel(index_dict=index_dict)
        Gamma = Gamma.convert_to_shorthand()

        g_inv = self.invert_metric(g=g, index_dict=index_dict)

        dim = len(index_dict)
        for i in range(dim):
            for j in range(dim):
                for k in range(dim):
                    r = index_dict[i] # rho
                    m = index_dict[j] # mu
                    n = index_dict[k] # nu  
                    rmn = r+m+n 
                    rs = sp.symbols(r)
                    ms = sp.symbols(m)
                    ns = sp.symbols(n)
                    for l in range(dim):
                        s = index_dict[l] 
                        rl = r+s
                        nl = n+s
                        lm = s+m
                        mn = m+n
                        ls = sp.symbols(index_dict[l])
                        try:
                            Gamma[rmn] += (1/2)*(g[rl])**(-1)*(sp.diff(g[nl], ms) + sp.diff(g[lm], ns) - sp.diff(g[mn], ls))
                            #Gamma[rmn] += (1/2)*g_inv[rl]*(sp.diff(g[nl], ms) + sp.diff(g[lm], ns) - sp.diff(g[mn], ls))
                        except ZeroDivisionError:
                            ppp = None
                            #print('g[{}] is zero, cannot evaluate.'.format(rl))
                    if simplify is True:
                        Gamma[rmn] = sp.simplify(Gamma[rmn])
        return Gamma

    def calculate_Riemann(self, Gamma, index_dict, simplify):
        """Calculate the Riemman tensor for given Christoffel symbols ``Gamma``. 

        Parameters: 
        -----------
        Gamma : dict, sympy.core.symbols.Symbol
            Construct the Christoffel symbols first by running self.calculate_Christoffel(). 
            Then pass to this function. 

        index_dict : dict, string 
            The names of the indices as strings. 

        simplify : bool 
            If true, simplify with sp.simplify. Else do nothing. 

        Returns: 
        -------- 
        riemann : dict, sympy.core.symbols.Symbol 
            A dictionary containing the Riemann tensor entries as sympy symbols. 
        """

        riemann = Riemann(index_dict=index_dict) 
        riemann = riemann.convert_to_shorthand()
        
        dim = len(index_dict)
        for i in range(dim):
            for j in range(dim):
                for k in range(dim):
                    for p in range(dim):
                        r = index_dict[i] # rho
                        s = index_dict[j] # sigma 
                        m = index_dict[k] # mu 
                        n = index_dict[p] # nu
                        rsmn = r+s+m+n 
                        rsn = r+s+n
                        rsm = r+s+m
                        sn = s+n
                        sm = s+m
                        ms = sp.symbols(m)
                        ns = sp.symbols(n)
                        riemann[rsmn] = sp.diff(Gamma[rsn], ms) - sp.diff(Gamma[rsm], ns) 
                        for l in range(dim):
                            l_str = index_dict[l] #lambda
                            lm = l_str+m 
                            ln = l_str+n 
                            rlm = r+l_str+m
                            rln = r+l_str+n
                            lsn = l_str+s+n
                            lsm = l_str+s+m
                            riemann[rsmn] += (Gamma[rlm]*Gamma[lsn] - Gamma[rln]*Gamma[lsm]) 
                            # Note: not immediately obvious that this follows the Einstein 
                            # convention when re-packed, but it does! 
                        if simplify is True:
                            riemann[rsmn] = sp.simplify(riemann[rsmn])
        return riemann
                            
    def calculate_Ricci_tensor(self, riemann, index_dict, simplify):
        """Calculate the Ricci tensor for a given Riemann tensor input. 

        Parameters:
        -----------
        riemann : dict, sympy.core.symbols.Symbols
            Pass this function the output of ``self.calculate_Riemann()``. 
        
        index_dict : dict, string 
            The name of the indices as strings 

        simplify : bool 
            If true, simplify with sp.simplify. Else do nothing. 

        Returns: 
        --------
        ricci : dict, sympy.core.symbols.Symbol 
        """

        ricci = Metric(index_dict=index_dict) # reuse this code, it is general for a (0,2) tensor
        ricci = ricci.convert_to_shorthand()

        dim = len(index_dict)
        for i in range(dim):
            for j in range(dim):
                m = index_dict[i] # mu 
                n = index_dict[j] # nu 
                mn = m+n 
                for l in range(dim):
                    l_str = index_dict[l] # lambda
                    lmln = l_str+m+l_str+n 
                    ricci[mn] += riemann[lmln] 

                if simplify is True:
                    ricci[mn] = sp.simplify(ricci[mn])
        return ricci 

    def calculate_Ricci_scalar(self, ricci_tensor, g, index_dict, simplify):
        """Calculate the Ricci scalar for a given Ricci tensor. 

        Parameters: 
        -----------

        ricci_tensor : dict, sympy.core.symbols.Symbol 
            Pass this function the output of ``self.calculate_Ricci_tensor()``. 

        g : dict, sympy.core.symbols.Symbol 
            Pass this function the metric to be studied. 

        index_dict : dict, string 
            A dictionary with the names of the indices as strings. 

        Returns: 
        --------
        ricci_scalar : sympy.core.symbols.Symbol 
        """

        ricci_scalar = 0.0 

        dim = len(index_dict) 
        for i in range(dim):
            for j in range(dim):
                m = index_dict[i] # mu 
                n = index_dict[j] # nu 
                mn = m+n
                try: 
                    ricci_scalar += (g[mn])**(-1)*ricci_tensor[mn]
                except ZeroDivisionError:
                    ppp = None
                    #print('g[{}] is zero, cannot evaluate.'.format(mn))
        if simplify is True:
            ricci_scalar = sp.simplify(ricci_scalar)
        #sp.cancel(ricci_scalar)
        ricci_scalar = {'':ricci_scalar} # convert to dictionary for handling - scalars have 0 indices  
        return ricci_scalar 

    def calculate_Einstein_tensor(self, ricci_tensor, ricci_scalar, g, index_dict, simplify):
        """Calculate the Einstein tensor given a Ricci tensor, Ricci scalar, and 
        the metric. 

        Parameters: 
        ----------- 
        ricci_tensor : dict, sympy.core.symbols.Symbol 
            The Ricci tensor to be studied. 

        ricci_scalar : sympy.core.symbols.Symbol
            The Ricci scalar, connected to the Ricci tensor passed in. 

        g : dict, sympy.core.symbols.Symbol 
            The metric to be studied. 

        index_dict : dict, string 
            The names of the indices, as strings. 

        Notes: 
        ------
        Pass this function the output of ``self.calculate_Ricci_tensor()``, 
        ``self.calculate_Ricci_scalar()``, and the metric. 

        Returns: 
        -------- 
        einstein : dict, sympy.core.symbols.Symbol 
        """

        einstein = Metric(index_dict=index_dict) # reuse code; general purpose for (0,2) tensor
        einstein = einstein.convert_to_shorthand()

        dim = len(index_dict)
        for i in range(dim):
            for j in range(dim):
                m = index_dict[i] # mu 
                n = index_dict[j] # nu 
                mn = m+n 
                einstein[mn] = ricci_tensor[mn] - (1/2)*ricci_scalar['']*g[mn] # scalars have 0 indices 
                if simplify is True:
                    einstein[mn] = sp.simplify(einstein[mn])
        return einstein

    def analyze_metric(self, g, index_dict, simplify):
        """Pass this function the metric to be analyzed. 
        It will compute: 
        (1) Christoffel symbols
        (2) Riemann tensor
        (3) Ricci tensor
        (4) Ricci scalar
        (5) Einstein tensor

        Parameters: 
        -----------
        g : dict, sympy.core.symbols.Symbol
            The metric, with components specified in the 
            shorthand way i.e. g['tt'] etc.

        index_dict : dict, string 
            The name of the indices as variables. 
            E.g. for Schwarzchild: {0:'t', 1:'r', 2:'theta', 3:'phi'}
            
        simplify : dict, bool 
            A dictionary of True or False booleans, with keys specified 
            by the name of the item to be simplified. 
            
        Returns: 
        --------
        Gamma : dict, sympy.core.symbols.Symbol
            The Christoffel symbols in the usual shorthand
            dict convention. 

        riemann : dict, sympy.core.symbols.Symbol
            The Riemann tensor components in the usual 
            shorthand dict convention. 

        ricci_tensor : dict, sympy.core.symbols.Symbol
            The Ricci tensor components in the usual 
            shorthand dict convention. 

        ricci_scalar : sympy.core.symbols.Symbol
            The Ricci scalar as a symbol. 

        einstein_tensor : dict, sympy.core.symbols.Symbol
            The Einstein tensor components in the usual 
            shorthand dict convention. 
        """

        print('Calculating Christoffel symbols...') 
        Gamma = self.calculate_Christoffel(g=g, index_dict=index_dict, simplify=simplify['Gamma'])
        print('Done.')
        print('Calculating Riemann tensor...') 
        riemann = self.calculate_Riemann(Gamma=Gamma, index_dict=index_dict, simplify=simplify['riemann'])
        print('Done.')
        print('Calculating Ricci tensor...')
        ricci_tensor = self.calculate_Ricci_tensor(riemann=riemann, index_dict=index_dict, simplify=simplify['ricci_tensor'])
        print('Done.')
        print('Calculating Ricci scalar...') 
        ricci_scalar = self.calculate_Ricci_scalar(ricci_tensor=ricci_tensor, g=g, index_dict=index_dict, simplify=simplify['ricci_scalar'])
        print('Done.') 
        print('Calculating Einstein tensor...') 
        einstein = self.calculate_Einstein_tensor(ricci_tensor=ricci_tensor, ricci_scalar=ricci_scalar, g=g, index_dict=index_dict, simplify=simplify['einstein'])
        print('Done.')

        return {'Gamma':Gamma, 'riemann':riemann, 'ricci_tensor':ricci_tensor,
                'ricci_scalar':ricci_scalar, 'einstein':einstein}

    def print_results(self, analysis_results):
        """Print out the results of ``self.analyze_metric()``. 
        
        Parameters:
        -----------
        analysis_resuts : dict, dict, sp.core.symbols.Symbol
            The output of ``self.analyze_metric()``. 
        """

        for key in analysis_results: 
            print('-----------------------------------------------------') 
            print('{key} ANALYSIS'.format(key=key))
            for index_string in analysis_results[key]: 
                print('\t{key}[{index_string}] = {expr}'.format(key=key,
                                                                index_string=index_string,
                                                                expr=sp.latex(analysis_results[key][index_string])))
            print('-----------------------------------------------------')

    def evaluate(self, name, symbol, point): 
        """Evaluate ``symbol`` at ``point`` with variables 
        specified by ``index_dict``. 

        Parameters:
        -----------
        name : string 
             A string identifying what to call the symbol. 

        symbol : dict, sympy.core.symbols.Symbol
             The symbol to be evaluated. Should have the shorthand structure
             even if it is a scalar. For example, to evaluate a metric at a 
             point, pass the shorthand g['tt'] etc. version. 

        point : dict, float 
             A dictionary for which the key is the variable as a 
             sympy.core.symbols.Symbol object and the value is the
             numerical value at which to evaluate that symbol. 

        Returns:
        --------
        evaluated_symbol : dict, float 
             A dictionary for which the key is whatever the key was in the 
             original ``symbol`` dict and the value is the value of that 
             ``symbol`` dict evaluated at ``point``. 
        """

        evaluated_symbol = symbol

        for key in evaluated_symbol:
            for var in point:
                #print('Evaluating {name}[{key}] at {var} = {point}.'.format(name=name,
                #key=key,
                #var=var,
                #point=point[var]))
                evaluated_symbol[key] = evaluated_symbol[key].subs(var, point[var])

        return evaluated_symbol

