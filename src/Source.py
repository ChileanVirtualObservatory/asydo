
class Source:
    """A source of EM Waves"""

    def __init__(self, log, name, alpha, delta):
        self.log = log
        log.write('Source \'' + name + '\' added \n')
        self.alpha = alpha
        self.delta = delta
        self.name = name
        self.comp = list()

    def addComponent(self, model):
        code = self.name + '-c' + str(len(self.comp) + 1) #+ '-r' + str(self.alpha) +'-d'+str(self.delta)
        self.comp.append(model)
        model.register(code,self.alpha,self.delta)
        self.log.write('Component added to \'' + self.name + '\' (' + code + ')\n')
        self.log.write('ModelInfo: '+ model.info() + '\n')
        
    def project(self, cube):
        for component in self.comp:
            self.log.write(' --> Component name: ' + component.comp_name + '\n')
            component.project(cube);
