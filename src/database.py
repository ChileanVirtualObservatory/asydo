import os
import sqlite3 as lite
import sys
import math
from astropy.io.votable.tree import Field as pField
from pprint import pprint



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

    def __init__(self, filename):
        self.name = filename

    def connect(self):
        try:
            self.pointer = lite.connect(self.name)
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

    def printTableDef(self,command):
        f = open("Table_def", "w")
        com, content = command.split("(")
        lines = content.split(",")
        print command
        print com
        map(pprint,lines)

    def genTable(self):
        command = "CREATE TABLE Catalog ("
        metadata = "CREATE TABLE Metadata (Column TEXT NOT NULL, Description TEXT NOT NULL)"
        insertMetadata = []
        count = 0
        for i in self.fields:
            name,description, dataType = i

            if count != 0:
                command = command + ", "

            command = command + " " + name.replace(" ", "_") + " " + dataType
            command2 = "INSERT INTO Metadata VALUES('" + name + "', '" + description + "')"
            insertMetadata.append(command2)
            count += 1

        command = command + ")"
        self.connect()

        self.pointer.execute(command)
        self.pointer.execute(metadata)
        for com in insertMetadata:
            self.pointer.execute(com)
        self.pointer.commit()

        self.Disconnect()

    def genInsertDataCommand(self, data, firstId):
        rawData = data._data
        id = firstId
        insertData = []
        for line in rawData:
            c = False
            command = "INSERT INTO Catalog VALUES("
            for value in line:
                if c:
                    command = command + ", "
                c = True
                if type(value).__name__ == "str":
                    temp = value.replace("'","''")
                    command = command + "'" + temp + "'"

                elif "int" in type(value).__name__ or "float" in type(value).__name__:
                    if math.isnan(value):
                        command = command + "'NaN'"
                    else:
                        command = command + str(value)
                elif "bool" in type(value).__name__:
                    if value:
                        command = command + str(1)
                    else:
                        command = command + str(0)
                else:
                    print type(value).__name__
                    print "Data Type not supported. Data not inserted in database. Exiting now!"
                    sys.exit()
            command = command + ")"
            insertData.append(command)
        return insertData

    def insertData(self,data):
        commands = self.genInsertDataCommand(data,0)
        self.connect()
        for com in commands:
            print com
            self.pointer.execute(com)
        self.pointer.commit()



    def deleteDB(self):
        os.remove(self.name)











