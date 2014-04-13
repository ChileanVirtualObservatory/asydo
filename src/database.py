import sqlite3 as lite
import sys
from astropy.io.votable.tree import Field as pField

SqlEquivalent = {
    "char":     "TEXT",
    "double":   "DOUBLE",
    "int":      "INT",
    "boolean":  "BOOLEAN"
}


class DataBase:
    name = ""
    connected = False
    fields = {}
    pointer = None

    def Connect(self, name):
        self.name = name
        try:
            self.pointer = lite.connect(name)
            self.connected = True
        except lite.Error, e:
            print "Error %s:" % e.args[0]
            sys.exit(1)

    def Disconnect(self):
        if self.pointer:
            self.pointer.close()

    def loadFields(self, fields):
        #recieves the fields from the VOtable and adds them to the class
        for f in fields:
            if isinstance(f, pField):
                name = f.name
                description = f.description
                type = SqlEquivalent[f.datatype]
                self.fields[name] = (description,type)





