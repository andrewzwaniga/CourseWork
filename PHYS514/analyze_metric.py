"""Analyze a given metric. 
"""

from metric_handler import * 

analysis = MetricAnalysis() 

problem_2 = False  
problem_3 = True 
problem_4 = False

if problem_2: 
    g_schwarzchild = analysis.init_schwarzchild()
    schwarzchild_analysis = analysis.analyze_metric(g=g_schwarzchild, index_dict=analysis.SCHWARZCHILD)
    analysis.print_results(g_schwarzchild, analysis.SCHWARZCHILD)

    g_deSitter = analysis.init_deSitter()
    deSitter_analysis = analysis.analyze_metric(g=g_deSitter, index_dict=analysis.DESITTER)
    analysis.print_results(g_deSitter, analysis.DESITTER)

    einstein_deSitter = deSitter_analysis['einstein'] 
    product = 0.0 
    for i in range(len(analysis.DESITTER)):
        for j in range(len(analysis.DESITTER)):
            m = analysis.DESITTER[i] 
            n = analysis.DESITTER[j] 
            mn = m+n
            try: 
                product += (g_deSitter[mn])**(-1)*einstein_deSitter[mn] 
            except ZeroDivisionError: 
                print('g[{mn}] is zero, cannot evaluate.'.format(mn=m+n))
    print('For de Sitter space, the product (g^mn)*(G_mn) is {product}.'.format(product=product))

if problem_3: 
    g_kerr = analysis.init_kerr()
    simplify = {'Gamma':True, 'riemann':False, 'ricci_tensor':False, 'ricci_scalar':False, 'einstein':False}
    kerr_analysis = analysis.analyze_metric(g=g_kerr, index_dict=analysis.KERR, simplify=simplify)
    analysis.print_results(kerr_analysis)
    


