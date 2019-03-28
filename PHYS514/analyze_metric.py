"""Analyze a given metric. 
"""

#from metric_handler import * 
from metric import Metric
from christoffel import Christoffel
from riemann import Riemann
from ricci import Ricci
from scalar import Scalar
from einstein import Einstein
from metric_analysis import MetricAnalysis

import sympy as sp 
import math 

analysis = MetricAnalysis() 

problem_2 = False
problem_3 = True
problem_4 = False

if problem_2: 
    print('#############################################################################################')
    print('##################################### Problem 2 #############################################')
    print('#############################################################################################')
    simplify = {'Gamma':True, 'riemann':True, 'ricci_tensor':True, 'ricci_scalar':True, 'einstein':True}
    g_schwarzchild = analysis.init_schwarzchild()
    schwarzchild_analysis = analysis.analyze_metric(g=g_schwarzchild, simplify=simplify)
    analysis.print_results(schwarzchild_analysis)

    g_deSitter = analysis.init_deSitter()
    deSitter_analysis = analysis.analyze_metric(g=g_deSitter, simplify=simplify)
    analysis.print_results(deSitter_analysis)

    einstein_deSitter = deSitter_analysis['einstein'] 
    g_deSitter_inv = analysis.invert_metric(g=g_deSitter)
    product = 0.0 
    for i in range(len(g_deSitter.index_dict)):
        for j in range(len(g_deSitter.index_dict)):
            m = g_deSitter.index_dict[i]  
            n = g_deSitter.index_dict[j] 
            mn = m+n
            try: 
                product += g_deSitter_inv.elements[mn]*einstein_deSitter.elements[mn] 
            except ZeroDivisionError: 
                ppp = 0
    print('For de Sitter space, the product (g^mn)*(G_mn) is {product}.'.format(product=product))

if problem_3:
    print('#############################################################################################')
    print('##################################### Problem 3 #############################################')
    print('#############################################################################################')
    
    latex_text_file = '/Users/azwaniga/Scratch/CourseWork/PHYS514/problem3_latex.txt'
    python_text_file = '/Users/azwaniga/Scratch/CourseWork/PHYS514/problem3_python.txt'
    t = sp.symbols('t')
    r = sp.symbols('r')
    theta = sp.symbols('theta')
    phi = sp.symbols('phi')
    a = sp.symbols('a')
    J = sp.symbols('J')
    M = sp.symbols('M')
    G = sp.symbols('G')

    g_kerr = analysis.init_kerr()
    g_kerr_inv = analysis.read_in_symbol_from_file(file_with_path=python_text_file, complex_symbol=g_kerr, symbol_type='metric')
    #g_kerr_inv = analysis.invert_metric(g_kerr)
    #styles = ['latex', 'python']
    #text_files = [latex_text_file, python_text_file]
    #for i in range(len(styles)):
    #    analysis.write_complex_symbol(file_with_path=text_files[i],
    #                                  complex_symbol=g_kerr_inv.elements,
    #                                  write_style=styles[i])
    simplify = {'Gamma':False, 'riemann':False, 'ricci_tensor':False, 'ricci_scalar':False, 'einstein':False}
    kerr_analysis = analysis.analyze_metric(g=g_kerr, g_inv=g_kerr_inv, simplify=simplify)
    einstein = kerr_analysis['einstein'] 
    point = {t:0, r:0.5, theta:math.pi/3, phi:math.pi/4, a:100, M:0.1, G:1} 
    print('Evaluating the Einstein tensor at this point:')
    for var in point:
        print('\t{var} = {point}'.format(var=var, point=point[var]))
    einstein_eval = analysis.evaluate(name='einstein', symbol=einstein, point=point)
    print('The result is...')
    for key in einstein_eval.elements:
        print('\tG[{key}] = {expr}'.format(key=key, expr=einstein_eval.elements[key]))
