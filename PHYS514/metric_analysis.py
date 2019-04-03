from metric import Metric
from christoffel import Christoffel
from riemann import Riemann
from ricci import Ricci
from scalar import Scalar
from einstein import Einstein
from ricci import Ricci

import sympy as sp
import numpy as np
import math 
    
class MetricAnalysis(object):
    """Class with functions that will operate on a given metric
    to produce the Christoffel symbols, the Riemann tensor, the Ricci
    tensor, Ricci scalar, and Einstein tensor. 
    """
    
    def __init__(self):
        self.SCHWARZCHILD = {0:'t', 1:'r', 2:'theta', 3:'phi'}
        self.S2 = {0:'theta', 1:'phi'}
        self.DESITTER = {0:'t', 1:'x', 2:'y', 3:'z'} 
        self.KERR = {0:'t', 1:'r', 2:'theta', 3:'phi'}
        self.FRW = {0:'t', 1:'r', 2:'theta', 3:'phi'}

    def init_schwarzchild(self):
        """Initialize Schwarzchild metric.
        """
        t = sp.symbols(self.SCHWARZCHILD[0]) 
        r = sp.symbols(self.SCHWARZCHILD[1])
        theta = sp.symbols(self.SCHWARZCHILD[2])
        phi = sp.symbols(self.SCHWARZCHILD[3])
        R = sp.symbols('R')

        g = Metric(index_dict=self.SCHWARZCHILD)
        g.convert_to_shorthand()
        g.elements['tt'] = 1 - R/r
        g.elements['rr'] = (1 - R/r)**(-1)
        g.elements['thetatheta'] = r**2
        g.elements['phiphi'] = r**2*(sp.sin(theta))**2

        return g 

    def init_schwarzchild_func(self):
        """Initialize Schwarzchild metric
        but the variables are now functions 
        of a parameter.
        """
        t = sp.Function(self.SCHWARZCHILD[0])
        r = sp.Function(self.SCHWARZCHILD[1])
        theta = sp.Function(self.SCHWARZCHILD[2])
        phi = sp.Function(self.SCHWARZCHILD[3])
        R = sp.symbols('R')
        y = sp.symbols('y')

        g = Metric(index_dict=self.SCHWARZCHILD)
        g.convert_to_shorthand()
        g.elements['tt'] = 1 - R*r(y)**(-1)
        g.elements['rr'] = (1-R*r(y)**(-1))**(-1)
        g.elements['thetatheta'] = r(y)**2
        g.elements['phiphi'] = r(y)**2*(sp.sin(theta(y)))**2

        return g 

    def init_S2(self):
        """Initialize the metric on the 2-sphere, S2. 
        """
        theta = sp.symbols(self.S2[0])
        phi = sp.symbols(self.S2[1])
        g = Metric(index_dict=self.S2)
        g.convert_to_shorthand()
        g.elements['thetatheta'] = 1 
        g.elements['phiphi'] = (sp.sin(theta))**2

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
        g.convert_to_shorthand()
        g.elements['tt'] = -1 
        g.elements['xx'] = sp.exp(2*H*t)
        g.elements['yy'] = g.elements['xx'] 
        g.elements['zz'] = g.elements['xx'] 

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
        g.convert_to_shorthand()

        g.elements['tt'] = -(1-2*G*M*r/rho**2)
        g.elements['tt'] = g.elements['tt'].subs(rho**2, rho2_expr)
        
        g.elements['tphi'] = -2*G*M*a*r*(sp.sin(theta))**2/rho**2
        g.elements['tphi'] = g.elements['tphi'].subs(rho**2, rho2_expr)
        g.elements['phit'] = g.elements['tphi'] 
        
        g.elements['rr'] = rho**2/Delta
        g.elements['rr'] = g.elements['rr'].subs(rho**2, rho2_expr)
        g.elements['rr'] = g.elements['rr'].subs(Delta, Delta_expr)
        
        g.elements['thetatheta'] = rho**2
        g.elements['thetatheta'] = g.elements['thetatheta'].subs(rho**2, rho2_expr)
                
        g.elements['phiphi'] = ((sp.sin(theta))**2/rho**2)*((r**2 + a**2)**2 - a**2*Delta*(sp.sin(theta))**2)
        g.elements['phiphi'] = g.elements['phiphi'].subs(rho**2, rho2_expr)
        g.elements['phiphi'] = g.elements['phiphi'].subs(Delta, Delta_expr)
        
        return g 

    def init_FRW(self): 
        """Initialize the FRW metric on 4D spacetime.
        """ 
        
        t = sp.symbols('t') 
        r = sp.symbols('r') 
        kappa = sp.symbols('k') 
        theta = sp.symbols('theta')
        phi = sp.symbols('phi') 
        a = sp.Function('a')

        g = Metric(index_dict=self.FRW) 
        g.convert_to_shorthand() 

        g.elements['tt'] = -1 
        g.elements['rr'] = a(t)**2/(1-kappa*r**2)
        g.elements['thetatheta'] = a(t)**2*r**2
        g.elements['phiphi'] = a(t)**2*r**2*sp.sin(theta)**2

        return g

    def invert_metric(self, g): 
        """ Invert a given metric g and return the inverted metric 
        as a shorthand dictionary of sympy.core.symbols.Symbol objects.
        """

        N = int(np.sqrt(len(g.elements))) 
        g_matrix = sp.MatrixSymbol('g_matrix', N, N)
        g_matrix = sp.Matrix(g_matrix)
        print('Converting the metric to sp.Matrix object...') 
        for i in range(N): 
            for j in range(N):
                m = g.index_dict[i]
                n = g.index_dict[j] 
                mn = m+n
                g_matrix[i,j] = g.elements[mn] 

        print('Inverting the matrix...') 
        g_matrix_inv = g_matrix**(-1)        
        g_inv = Metric(index_dict=g.index_dict)
        for i in range(N):
            for j in range(N): 
                m = g_inv.index_dict[i]
                n = g_inv.index_dict[j] 
                mn = m+n
                g_inv.elements[mn] = g_matrix_inv[i, j]  
        print('Succesfully inverted the matrix.') 
        return g_inv

    def calculate_Christoffel(self, g=None, g_inv=None, simplify=None): 
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

        Gamma = Christoffel(index_dict=g.index_dict)
        Gamma.convert_to_shorthand()
        
        if g_inv == None:
            print('Trying to invert the metric...') 
            g_inv = self.invert_metric(g=g)
            print('Finished inverting the metric.')

        dim = len(Gamma.index_dict)
        for i in range(dim):
            for j in range(dim):
                for k in range(dim):
                    r = Gamma.index_dict[i] # rho
                    m = Gamma.index_dict[j] # mu
                    n = Gamma.index_dict[k] # nu  
                    rmn = r+m+n 
                    rs = sp.symbols(r)
                    ms = sp.symbols(m)
                    ns = sp.symbols(n)
                    for l in range(dim):
                        s = Gamma.index_dict[l] 
                        rl = r+s
                        nl = n+s
                        lm = s+m
                        mn = m+n
                        ls = sp.symbols(Gamma.index_dict[l])
                        try:
                            Gamma.elements[rmn] += (1/2)*g_inv.elements[rl]*(sp.diff(g.elements[nl], ms) + sp.diff(g.elements[lm], ns) - sp.diff(g.elements[mn], ls))
                        except ZeroDivisionError:
                            ppp = None
                        if simplify is True:
                            Gamma.elements[rmn] = sp.simplify(Gamma.elements[rmn])
        return Gamma

    def calculate_Riemann(self, Gamma, simplify):
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

        riemann = Riemann(index_dict=Gamma.index_dict) 
        riemann.convert_to_shorthand()
        
        dim = len(riemann.index_dict)
        for i in range(dim):
            for j in range(dim):
                for k in range(dim):
                    for p in range(dim):
                        r = riemann.index_dict[i] # rho
                        s = riemann.index_dict[j] # sigma 
                        m = riemann.index_dict[k] # mu 
                        n = riemann.index_dict[p] # nu
                        rsmn = r+s+m+n 
                        rsn = r+s+n
                        rsm = r+s+m
                        sn = s+n
                        sm = s+m
                        ms = sp.symbols(m)
                        ns = sp.symbols(n)
                        riemann.elements[rsmn] = sp.diff(Gamma.elements[rsn], ms) - sp.diff(Gamma.elements[rsm], ns) 
                        for l in range(dim):
                            l_str = riemann.index_dict[l] #lambda
                            lm = l_str+m 
                            ln = l_str+n 
                            rlm = r+l_str+m
                            rln = r+l_str+n
                            lsn = l_str+s+n
                            lsm = l_str+s+m
                            riemann.elements[rsmn] += (Gamma.elements[rlm]*Gamma.elements[lsn] - Gamma.elements[rln]*Gamma.elements[lsm]) 
                            # Note: not immediately obvious that this follows the Einstein 
                            # convention when re-packed, but it does! 
                        if simplify is True:
                            riemann.elements[rsmn] = sp.simplify(riemann.elements[rsmn])
        return riemann
                            
    def calculate_Ricci_tensor(self, riemann, simplify):
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

        ricci = Ricci(index_dict=riemann.index_dict) 
        ricci.convert_to_shorthand()

        dim = len(ricci.index_dict)
        for i in range(dim):
            for j in range(dim):
                m = ricci.index_dict[i] # mu 
                n = ricci.index_dict[j] # nu 
                mn = m+n 
                for l in range(dim):
                    l_str = ricci.index_dict[l] # lambda
                    lmln = l_str+m+l_str+n 
                    ricci.elements[mn] += riemann.elements[lmln] 

                if simplify is True:
                    ricci.elements[mn] = sp.simplify(ricci.elements[mn])
        return ricci 

    def calculate_Ricci_scalar(self, ricci_tensor=None, g=None, g_inv=None, simplify=None):
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
        if g_inv == None:
            print('Trying to invert the metric...')
            g_inv = self.invert_metric(g=g)
            print('Succesfully inverted the metric.') 

        ricci_scalar = Scalar(index_dict=ricci_tensor.index_dict)

        dim = len(ricci_tensor.index_dict) 
        for i in range(dim):
            for j in range(dim):
                m = ricci_tensor.index_dict[i] # mu 
                n = ricci_tensor.index_dict[j] # nu 
                mn = m+n
                try: 
                    ricci_scalar.elements[''] += g_inv.elements[mn]*ricci_tensor.elements[mn]
                except ZeroDivisionError:
                    ppp = None
            if simplify is True:
                ricci_scalar.elements[''] = sp.simplify(ricci_scalar.elements[''])
        return ricci_scalar 

    def calculate_Einstein_tensor(self, ricci_tensor, ricci_scalar, g, simplify):
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

        einstein = Einstein(index_dict=ricci_tensor.index_dict) 
        einstein.convert_to_shorthand()

        dim = len(einstein.index_dict)
        for i in range(dim):
            for j in range(dim):
                m = einstein.index_dict[i] # mu 
                n = einstein.index_dict[j] # nu 
                mn = m+n 
                einstein.elements[mn] = ricci_tensor.elements[mn] - (1/2)*ricci_scalar.elements['']*g.elements[mn] # scalars have 0 indices 
                if simplify is True:
                    einstein.elements[mn] = sp.simplify(einstein.elements[mn])
        return einstein

    def solve_geodesic(self, g, likeness):
        """Calculate and solve the geodesic equations
        for a static, spherically symmetric metric for a timelike, 
        null, or spacelike worldline that is parameterized 
        by an affine parameter, in order to arrive at an 
        expression for the effective potential in terms of a 
        radial coordinate.  

        Parameters: 
        -----------
        g : Metric 
            Note that this function should be passed a Metric instance
            but the entries of g.elements should be sp.Function objects
            and not sp.symbols objects, otherwise the variables in the 
            geodesic will not coincide with those in the metric! 

        likeness : string 
            One of 'timelike', 'spacelike', or 'null' in order to 
            determine a value for ``epsilon`` that will enter into
            the RHS of the equation 
            \[
            \epsilon = -g_{\mu\nu}\frac{dx^\mu}{d\lambda}\frac{dx^\nu}{d\lambda}
            ]]
        """

        worldline_likeness = {'timelike':-1, 'spacelike':1, 'null':1}
        epsilon = worldline_likeness[likeness]  

        g_inv = self.invert_metric(g=g) 
        Gamma = self.calculate_Christoffel(g=g, g_inv=g_inv, simplify=True) 
        
        t = sp.Function(g.index_dict[0])
        r = sp.Function(g.index_dict[1])
        theta = sp.Function(g.index_dict[2])
        phi = sp.Function(g.index_dict[3])
        
        y = sp.symbols('y')
        
        worldline = {'t':t, 'r':r, 'theta':theta, 'phi':phi}

        K = len(Gamma.index_dict)
        sum = 0 # for geodesic equation
        ds2 = 0 # for metric equation $ds^2 = g_{\mu\nu}\dot{x}^\mu\dot{x}^\nu = -1$ for timelike
        for i in range(K):
            for j in range(K):
                m = Gamma.index_dict[i] # mu
                n = Gamma.index_dict[j] # nu 
                rmn = 'r'+m+n
                mn = m+n
                xm = worldline[m] 
                xn = worldline[n] 
                sum += Gamma.elements[rmn]*(xm(y).diff(y))*(xn(y).diff(y))
                ds2 += g.elements[mn]*(xm(y).diff(y))*(xn(y).diff(y))
        geodesic = sp.Eq(r(y).diff(y,y)+sum, 0)
        ds2 = sp.Eq(ds2, epsilon)  

        # Now fix theta=pi/2 because the direction of 
        # angular momentum is conserved (there are two 
        # Killing vectors d_t and d_phi for a static,
        # spherically symmetric metric) 

        geodesic = geodesic.subs(theta(y).diff(y), 0)
        geodesic = geodesic.subs(theta(y), math.pi/2)
        ds2 = ds2.subs(theta(y).diff(y), 0)
        ds2 = ds2.subs(theta(y), math.pi/2) 
        ds2 = ds2.evalf()
        print(sp.latex(geodesic))
        print(sp.latex(ds2))

        geodesic = sp.dsolve(geodesic)
        L = sp.symbols('L')
        E = sp.symbols('E')
        R = sp.symbols('R')
        ds2 = ds2.subs(phi(y).diff(y), L*r(y)**(-2))
        ds2 = ds2.subs(t(y).diff(y), E*(1-R*r(y)**(-1))**(-1))
        ds2 = sp.solve(ds2, (r(y).diff(y))**2)
        print('The solution for $\Big{(}\\frac{dr}{dy}\Big{)}^2$ is:')
        print(sp.latex(ds2))
        # At this point there should be just ddot(r) and ddot(t) 
        # left between these two equations and we should be 
        # able to sub ddot(r) into ddot(t) equation and then 
        # integrate to get dot(t) and plug it back into ddot(r) etc. 


    def analyze_metric(self, g, g_inv, simplify):
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
        Gamma = self.calculate_Christoffel(g=g, g_inv=g_inv, simplify=simplify['Gamma'])
        print('Done.')
        print('Calculating Riemann tensor...') 
        riemann = self.calculate_Riemann(Gamma=Gamma, simplify=simplify['R'])
        print('Done.')
        print('Calculating Ricci tensor...')
        ricci_tensor = self.calculate_Ricci_tensor(riemann=riemann, simplify=simplify['P'])
        print('Done.')
        print('Calculating Ricci scalar...') 
        ricci_scalar = self.calculate_Ricci_scalar(ricci_tensor=ricci_tensor, g=g, g_inv=g_inv, simplify=simplify['Q'])
        print('Done.') 
        print('Calculating Einstein tensor...') 
        einstein = self.calculate_Einstein_tensor(ricci_tensor=ricci_tensor, ricci_scalar=ricci_scalar, g=g, simplify=simplify['G'])
        print('Done.')

        return {'Gamma':Gamma, 'R':riemann, 'P':ricci_tensor,
                'Q':ricci_scalar, 'G':einstein}

    def print_results(self, analysis_results):
        """Print out the results of ``self.analyze_metric()``. 
        
        Parameters:
        -----------
        analysis_resuts : dict, dict, sp.core.symbols.Symbol
            The output of ``self.analyze_metric()``. 
        """

        for key in analysis_results: 
            print('-----------------------------------------------------')
            print('SHOWING ONLY NONZERO COMPONENTS...')
            print('{key} ANALYSIS'.format(key=key))
            for index_string in analysis_results[key].elements: 
                if analysis_results[key].elements[index_string] != 0:
                    print('\t{key}_{index_string} &= {expr}\\\\'.format(key=key,
                                                                    index_string=index_string,
                                                                    expr=sp.latex(analysis_results[key].elements[index_string])))
            print('ALL OTHER COMPONENTS ARE ZERO.')
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

        for key in evaluated_symbol.elements:
            for var in point:
                evaluated_symbol.elements[key] = evaluated_symbol.elements[key].subs(var, point[var])

        return evaluated_symbol

    def write(self, file_with_path, data):
        """Write important calculations to a file.
        Be sure to hand this function string formatted `data` 
        that is suitable for writing to a text file. 

        Parameters:
        -----------
        file_with_path : string 
            The file name with full path specified, to which
            the information will be written. 

        data : string
            The contents to write to the file. Note that only
            string type is supported because we only support 
            writing to a text file currently. 

        Returns:
        --------
        """

        text_file = open(file_with_path, 'a+')
        text_file.write('\n')
        text_file.write('\t')
        text_file.write(data)
        text_file.write('\n')
        text_file.close()

    def write_symbol(self, file_with_path, symbol, write_style):
        """Write either a LaTeX or Python string for ``symbol`` 
        to a text file.
        
        Parameters:
        -----------
        file_with_path : string 
            The file name with full path specified, to which 
            the information will be written. 

        symbol : sympy.core.symbols.Symbol 
            The symbol to be written. 

        write_style : string 
            If you want to write a LaTeX string or a Python-compatible 
            string (meaning the string could be copy-pasted into Python
            as the value of a sympy expression) you can specify this 
            by passing 'latex' or 'python'. 

        Returns:
        --------
        """

        data = None 
        if write_style == 'latex': 
            data = sp.latex(symbol)
        elif write_style == 'python': 
            data = str(symbol) 
        else: 
            raise TypeError("Error writing symbol: specify either 'latex' or 'python' for `write_style`.")
 
        self.write(file_with_path=file_with_path, data=data)

    def write_complex_symbol(self, file_with_path, complex_symbol, write_style):
        """Write a complex symbol to a text file as either a LaTeX or Python string.
        A complex symbol means that a dictionary is expected to be passed to 
        this function, where the keys are labels for the symbol that is contained 
        in the corresponding value for the key. 

        The purpose is to accomodate writing the dictionary-styled symbols that are 
        used in this analysis script. Namely, the ``elements`` attribute of the classes 
        called ``Metric``, ``Christoffel``, ``Riemann``, ``Ricci``, ``Scalar``, and 
        ``Einstein``. 
        """
        
        for key in complex_symbol:
            with open(file_with_path, 'a+') as f:
                s = '\t[{key}]'.format(key=key)
                f.write(s)
            self.write_symbol(file_with_path=file_with_path, 
                              symbol=complex_symbol[key],
                              write_style=write_style)
        
    def read_in_symbol_from_file(self, file_with_path, complex_symbol, symbol_type):
        from sympy.parsing.sympy_parser import parse_expr 
        """Read a text file that contains a string version 
        of a complex symbol and return a complex symbol 
        with the entries that were stored in the text file. 
        Assumes that the text file has been created with 
        ``self.write_complex_symbol()``; i.e. this function
        is not fit for reading arbitrary text files. 
        """

        index_list = []
        if symbol_type == 'metric':
            constructed_symbol = Metric(index_dict=complex_symbol.index_dict)
        elif symbol_type == 'christoffel':
            constructed_symbol = Christoffel(index_dict=complex_symbol.index_dict)
        elif symbol_type == 'riemann':
            constructed_symbol = Riemann(index_dict=complex_symbol.index_dict)
        elif symbol_type == 'ricci':
            constructed_symbol = Ricci(index_dict=complex_symbol.index_dict)
        elif symbol_type == 'einstein':
            constructed_symbol = Einstein(index_dict=complex_symbol.index_dict)
        else:
            raise TypeError('Got a bad symbol type when trying to read in a symbol from a text file.') 
        
        constructed_symbol.convert_to_shorthand()

        with open(file_with_path, 'r') as text_file:
            text = text_file.read() 
            for key in complex_symbol.elements:
                index_list.append(key) 
            l = len(index_list)
            #print('Here is the index list: {list}.'.format(list=index_list))
            for i in range(l):
                s = '[{ss}]'.format(ss=index_list[i])
                #print('Searching for the end of {s}...'.format(s=s))
                start = text.find(s) + len(s) # start of the symbolic string 
                #print('Found the end of {s} at {loc}.'.format(s=s, loc=start))
                if i + 1 < l:
                    t = '[{tt}]'.format(tt=index_list[i+1])
                    #print('Searching for {t}...'.format(t=t))
                    end = text.find(t) 
                    #print('Found {t} at {loc}.'.format(t=t, loc=end))
                    symbol_string = text[start:end]
                    #print('Here is the symbol we found:\n')
                    #print('\t' + symbol_string)
                    output_symbol = parse_expr(symbol_string, evaluate=False)
                else:     
                    # handle the last symbol, which won't be followed 
                    # by another ['xx'] flag...
                    symbol_string = text[start:] # go to end of the string!
                    #print('Here is the symbol we found:\n')
                    #print('\t' + symbol_string)
                    output_symbol = parse_expr(symbol_string, evaluate=False)
                #print('Adding the following symbol to [{key}]:'.format(key=index_list[i]))
                constructed_symbol.elements[index_list[i]] = output_symbol 
                #print(constructed_symbol.elements[index_list[i]])

            return constructed_symbol
        
