from Source import *
from Cube import *


class Universe:
    """A synthetic universe where to put synthetic objects."""

    def __init__(self, log):
        self.log=log
        self.sources = dict()

    def createSource(self, name, alpha , delta):
        self.sources[name] = Source(self.log, name, alpha, delta)

    def addComponent(self, source_name, model):
        self.sources[source_name].addComponent(model)

    def genCube(self, name, cubespec, filename):
        cube = Cube(self.log, name, cubespec)
        for src in self.sources: 
            self.log.write('   * Source: ' + src + '\n')
            self.sources[src].project(cube)
        self.log.write('   * Saving FITS: ' + filename + '\n')
        cube.saveFits(self.sources, filename)
        return cube

    def removeSource(self, name):
        self.log.write('Removing source '+name)
        return self.sources.remove(name)

