"""
Generate an encyclopedia of numerical methods, using NodePy.
The encyclopedia may include:

    - A page for each method, with its coefficients and properties listed and plotted
    - A system of automatically generated tags, with an index page of tags and a page for each tag listing all so-tagged methods
    - More advanced: the ability to generate method-to-method comparisons
    - The ability to search for all methods satisfying a given set of criteria
"""
from nodepy import rk, lm, lsrk
import os

parent_path = os.path.abspath(os.path.dirname(os.path.realpath(__file__))+'/..')
print parent_path
img_dir = 'img'
methods_dir = 'methods'
top_dir='numipedia'
top_dir_abs = '/'+top_dir

def write_numipedia(rewrite_method_pages=True):
    """
    Main function to write the whole numipedia site.
    Loops over all methods in rk.loadRKM('All') and writes a page for each.

    TODO:

        - Add parameterized method families
        - Add more methods
    """
    methods = rk.loadRKM('All')

    for k in range(2,6):
        method = lm.Adams_Bashforth(k)
        methods[method.name] = method
        method = lm.Adams_Moulton(k)
        methods[method.name] = method
        method = lm.backward_difference_formula(k)
        methods[method.name] = method

    for key in ['DDAS47','LDDC46','RK45[2R]C','RK58[3R]C','RK59[2R]C']:
        method = lsrk.load_2R(key)
        methods[method.name] = method

    for key in ['NSSP32','NSSP33']:
        del methods[key]

    method = rk.SSPRK2(4)
    methods[method.name] = method
    method = rk.SSPRK3(4)
    methods[method.name] = method
    method = rk.SSPRK3(9)
    methods[method.name] = method
    method = rk.SSPIRK3(3)
    methods[method.name] = method
    method = rk.SSPIRK3(4)
    methods[method.name] = method

    method = rk.RKC1(5)
    methods[method.name] = method
    method = rk.RKC1(20)
    methods[method.name] = method
    method = rk.RKC2(5)
    methods[method.name] = method
    method = rk.RKC2(20)
    methods[method.name] = method

    method = rk.DC_pair(3)
    methods[method.name] = method
    method = rk.DC_pair(4)
    methods[method.name] = method

    method = rk.extrap_pair(4,base='euler')
    methods[method.name] = method
    method = rk.extrap_pair(6,base='euler')
    methods[method.name] = method

    method = rk.extrap_pair(2,base='midpoint')
    methods[method.name] = method
    method = rk.extrap_pair(3,base='midpoint')
    methods[method.name] = method
    method = rk.extrap_pair(4,base='midpoint')
    methods[method.name] = method

    write_index_page(methods,index_file_path=os.path.join(parent_path,top_dir))

    if rewrite_method_pages:
        methods_path = os.path.join(parent_path,top_dir,methods_dir)
        print methods_path
        if not os.path.exists(methods_path):
            os.makedirs(methods_path)
        img_path = os.path.join(parent_path,top_dir,img_dir)
        print img_path
        if not os.path.exists(img_path):
            os.makedirs(img_path)

        template_file = os.path.join(parent_path,top_dir,'method_template.html')
        for key,method in methods.iteritems():
            fname = os.path.join(methods_path,method.shortname+'.html')
            print fname
            s = method_page(method,template_file)

            outfile = open(fname,'w')
            outfile.write(s)
            outfile.close()


def method_page(method,template_file='method_template.html'):
    """
    Writes an HTML page for a specified method, using the template file.
    Uses mako to do simple substitutions.

    TODO:

        - Add more method properties
        - Improve page layout (in template)
    """
    from mako.template import Template
    import matplotlib.pyplot as plt

    # Plot stability region.
    plot_file_name = method.shortname+'.png'
    plot_file = os.path.join(parent_path,top_dir,img_dir,plot_file_name)
    plot_urlpath = os.path.join(top_dir_abs,img_dir,plot_file_name)
    fig = method.plot_stability_region(to_file=plot_file,longtitle=False)
    plt.close()

    if isinstance(method,rk.RungeKuttaMethod):
        # Compute and render stability function
        p,q=method.stability_function()
        from sympy import symbols, latex, factorial, Rational

        # Make order-correct coeffs exact
        for pol in (p,q):
            for i, coeff in enumerate(pol.c[::-1]):
                if abs(coeff-1./factorial(i))<1.e-14:
                    pol.c[-i-1] = Rational(1,factorial(i))

        z = symbols('z')
        pp = sum(co*z**i for i,co in enumerate(p.c[::-1]))
        qq = sum(co*z**i for i,co in enumerate(q.c[::-1]))
        stabfun = "$$"+latex(pp/qq,order='old')+"$$"
        stage_order = tex(method.__num__().stage_order())
    else:
        stabfun = None
        stage_order = method.order()

    if method.info is '' and hasattr(method,'mtype'):
        method.info = method.mtype

    mytemplate = Template(filename=template_file)
    return mytemplate.render(method=method,
                          name=method.name,
                          desc=method.info,
                          plot_urlpath=plot_urlpath,
                          butcher=method.latex(),
                          stabfun=stabfun,
                          amrad=tex(method.absolute_monotonicity_radius()),
                          order = tex(method.p),#__num__().order()),
                          stage_order = stage_order
                          )


def write_index_page(methods,fname='index.html',
                     template_file='index_template.html',
                     index_file_path='.',
                     methods_urlpath='/numipedia/methods'):
    from mako.template import Template
    mytemplate = Template(filename = os.path.join(index_file_path,template_file))

    method_props = {}
    for method_name, method in methods.iteritems():
        method_props[method.shortname] = {}
        methdict = method_props[method.shortname]
        properties = ["mix"]

        if method.is_explicit(): properties.append("explicit")
        else: properties.append("implicit")

        if isinstance(method,lm.LinearMultistepMethod):
            properties.append('multistep')
        elif isinstance(method,rk.RungeKuttaMethod):
            properties.append('runge-kutta')
            if hasattr(method,'mtype'):
                if method.mtype == 'Diagonally implicit Runge-Kutta method':
                    properties.append('diagonally-implicit')
        if isinstance(method,lsrk.TwoRRungeKuttaMethod) or \
           isinstance(method,lsrk.TwoSRungeKuttaMethod) or \
           isinstance(method,lsrk.TwoSRungeKuttaPair):
            properties.append('low-storage')

        if method.absolute_monotonicity_radius()>1.e-10:
            properties.append("ssp")

        properties.append("order-%s" % str(method.p))

        if hasattr(method,'embedded_method'):
            properties.append("pair")
        else:
            properties.append("not-pair")

        methdict['class string'] = " ".join(properties)
        methdict['name'] = method.name

    s = mytemplate.render(method_props=method_props,methods_urlpath=methods_urlpath)
    with open(os.path.join(index_file_path,fname),'w') as outfile:
        outfile.write(s)

def tex(s):
    from nodepy.snp import printable
    return "$"+printable(s,return_zero=True)+"$"

if __name__ == "__main__":
    write_numipedia()
