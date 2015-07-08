from Jinja2 import Environment, PackageLoader

env = Environment(loader=Packageloader('/users/asood/My Documents', 'Templates'))

template = env.get_template('test_template.html')
