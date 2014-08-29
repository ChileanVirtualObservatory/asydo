
class Source:
    """A source of EM Waves"""

    def __init__(self, log, name, alpha, delta):
        """ Parameters:
               * log: logging descriptor
               * name: a name of the source
               * alpha: right ascencion 
               * delta: declination
        """
        self.log = log
        log.write('+++ Source \'' + name + '\' added\n')
        self.alpha = alpha
        self.delta = delta
        self.name = name
        self.comp = list()

    def addComponent(self, model):
        """ Defines a new component from a model
            Parameters:
               * model : a valid Component 
        """
        code = self.name + '-c' + str(len(self.comp) + 1) #+ '-r' + str(self.alpha) +'-d'+str(self.delta)
        self.comp.append(model)
        model.register(code,self.alpha,self.delta)
        self.log.write(' |- Component added to \'' + self.name + '\' (' + code + ')\n')
        self.log.write(' ---+ Model: '+ model.info() + '\n')
        
    def project(self, cube):
        """Projects all components in the Source"""
        for component in self.comp:
            self.log.write('  |- Projecting ' + component.comp_name + '\n')
            component.project(cube);
