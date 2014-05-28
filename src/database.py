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
    fields = []
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
            self.connected = False

    def loadFields(self, fields):
        #recieves the fields from the VOtable and adds them to the class
        for f in fields:
            if isinstance(f, pField):
                name = f.name
                description = f.description
                type = SqlEquivalent[f.datatype]
                self.fields.append((name,description,type))

    def genTable(self):
        command = "CREATE TABLE Catalogo (ID INT PRIMARY KEY NOT NULL,"
        metadata = "CREATE TABLE Metadata (ID INT PRIMARY KEY NOT NULL, Column TEXT NOT NULL, Description TEXT NOT NULL)"
        insertmetadata = []
        contador = 0
        for i in self.fields:
            name,descripcion, tipodato = i

            if contador != 0:
                command = command + ", "


            command = command + " " + name.replace(" ", "_") + " " + tipodato
            command2 = "INSERT INTO Metadata VALUES(" + str(contador) + ", '" + name + "', '" + descripcion + "')"
            insertmetadata.append(command2)
            contador += 1

        command = command + ")"
        self.Connect("ASYDOGet.DB")

        self.pointer.execute(command)
        self.pointer.execute(metadata)
        for com in insertmetadata:
            self.pointer.execute(com)
        self.pointer.commit()

        self.Disconnect()










