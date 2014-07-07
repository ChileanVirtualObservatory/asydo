#SPEED_OF_LIGHT = 299792458.0
#S_FACTOR = 2.3548200450309493820231386529193992754947713787716410 #sqrt(8*ln2) #TODO Preguntar Profe que es esto?


class Source:
    """A source of EM Waves"""

    def __init__(self, log, name, alpha, delta):
        log.write('Source \'' + name + '\' added \n')
        self.log = log
        self.alpha = alpha
        self.delta = delta
        self.name = name
        self.comp = list()

    def addComponent(self, model):
        code = self.name + '-c' + str(len(self.comp) + 1 + '-r' + str(alpha) +'-d'+str(delta))
        self.comp.append(model)
        mode.register(self,code)
        log.write('Component added to \'' + self.name + '\' (' + code + ')\n')
        log.write('ModelInfo: '+ model.info() + '\n')
        
    def project(self, cube):
        for component in self.comp:
            log.write(' --> Component name: ' + struct.code + '\n')
            component.project(cube);
