import os
import sqlite3 as lite
import sys
import math
from astropy.io.votable.tree import Field as pField
from astropy.io.votable import parse_single_table
from astropy.io.votable.tree import Table as pTable
import urllib2
import csv



SqlEquivalent = {
    "char":     "TEXT",
    "double":   "DOUBLE",
    "int":      "INT",
    "boolean":  "BOOLEAN"
}

#USAGE EXAMPLE for single table
#
# location = './votables/band2.xml'
# db = DataBase("ASYDO.sqlite")
# db.loadVoTable(location)
#
class DataBase:
    name = ""
    connected = False
    fields = []
    pointer = None

    slap_serv = 'https://find.nrao.edu/splata-slap/slap'
    alma_band_freq = {'3': [88, 116], '4': [125, 163], '6': [211, 275], '7': [275, 373], '8': [385, 500], '9': [602, 720]}

    def __init__(self, dbName):
        self.name = dbName

    def connect(self):
        try:
            self.pointer = lite.connect(self.name+".sqlite")
            self.connected = True
        except lite.Error, e:
            print "Error %s:" % e.args[0]
            sys.exit(1)

    def disconnect(self):
        if self.pointer:
            self.pointer.close()
            self.connected = False

    def VOGetLines(self,log, prefix):
        c = 299792458.0
        for band in self.alma_band_freq:
            log.write('band=' + band + '\n')
            w_init = c / (self.alma_band_freq[band][0] * 1000000000.0)
            w_end = c / (self.alma_band_freq[band][1] * 1000000000.0)
            data = '?REQUEST=queryData&WAVELENGTH=' + \
                str(w_init) + '/' + str(w_end) + '&VERB=3'
            curl = self.slap_serv + data.encode('utf-8')
            log.write('  -> Downloading lines via SLAP:\n')
            log.write('  -> ' + curl + '\n')
            req = urllib2.Request(curl)
            response = urllib2.urlopen(req)
            votable = response.read()
            stri = prefix + band
            location = './votables/' + stri + '.xml'
            f = open(location, 'w')
            f.write(votable)
            f.close()

    def loadFields(self, fields):
        #recieves the fields from the VOtable and adds them to the class
        for f in fields:
            if isinstance(f, pField):
                name = f.name
                description = f.description
                type = SqlEquivalent[f.datatype]
                self.fields.append((name,description,type))

    def printTableDef(self,command):
        #Prints a definition table for the catalog into a txt file
        f = open(self.name+"-Table_def", "w")
        columns = ["Name","Description","DataType"]
        f.write(columns[0].rjust(20) +"\t | \t"+columns[1].ljust(100," ")+ "\t | \t" +columns[2]+"\n")
        divisor = "".rjust(25,"_")+"|"+ "".rjust(111,"_")+"|"+ "".rjust(20,"_")+"\n"
        f.write(divisor)
        for line in command:
            title, description, dataType = line
            f.write(title.rjust(20," ") +"\t | \t"+ description.ljust(100," ") +"\t | \t"+ dataType+"\n")
        f.close()

    def genTable(self):
        #Generates the Catalog and Metada Tables, and populates the latter.
        command = "CREATE TABLE Catalog ("
        metadata = "CREATE TABLE Metadata (Column TEXT NOT NULL, Description TEXT NOT NULL)"
        insertMetadata = []
        output_command = []
        count = 0
        for i in self.fields:
            name,description, dataType = i

            if count != 0:
                command += ", "

            command = command + " " + name.replace(" ", "_") + " " + dataType
            command2 = "INSERT INTO Metadata VALUES('" + name + "', '" + description + "')"
            o_command = (name.replace(" ", "_"), description, dataType)


            insertMetadata.append(command2)
            output_command.append(o_command)
            count += 1

        command = command + ")"
        self.printTableDef(output_command)
        self.connect()

        self.pointer.execute(command)
        self.pointer.execute(metadata)
        for com in insertMetadata:
            self.pointer.execute(com)
        self.pointer.commit()

        self.disconnect()

    def genInsertDataCommand(self, data):
        #Generates the commands for SQL Insertion
        rawData = data._data
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
        #inserts the data into the database
        commands = self.genInsertDataCommand(data)
        self.connect()
        for com in commands:
            self.pointer.execute(com)
        self.pointer.commit()

    def loadVoTable(self,location):
        #Generates the tables in the DB and loads the data.
        tbl = parse_single_table(location)
        if isinstance(tbl,pTable):
        # tbl.array contiene los datos
        # tbl.field contiene la metadata
            self.loadFields(tbl.fields)
            self.genTable()
            self.insertData(tbl.array)

    def loadMultipleVoTables(self,locations):
        #Loads the data of multiple VOTables, generating the definitions from the first VOTable in the location list
        #It assumes every VOTables has the same columns.
        tbl = parse_single_table(locations[0])
        if isinstance(tbl,pTable):
        # tbl.array contiene los datos
        # tbl.field contiene la metadata
            self.loadFields(tbl.fields)
            self.genTable()
            for place in locations:
                currentTable = parse_single_table(place)
                self.insertData(currentTable.array)


    def createDBFromCSV(self, filename):
        self.Connect()
        create = "CREATE TABLE Lines(ID INT PRIMARY KEY NOT NULL,SPECIES TEXT,CHEM_NAME TEXT,FREQ REAL,INTENSITY REAL,EL REAL)"
        drop = "DROP TABLE Lines"
        with open(filename, 'rb') as csvfile:
            sreader = csv.reader(csvfile, delimiter=':', quotechar='|')
            counter=0
            for row in sreader:
                if counter == 0:
                    try:
                       self.pointer.execute(drop)
                    except lite.OperationalError:
                       print "WARNING: Drop failed\n"
                    self.pointer.execute(create)
                    counter+=1
                else:
                    if len(row) < 10:
                       continue
                    species=row[0].replace("'","-")
                    chname=row[1].replace("'","-")
                    freq=row[2]
                    if freq=='':
                        freq=row[4]
                    insert="INSERT INTO Lines VALUES("+str(counter)+",'"+species+"','"+chname+\
                           "',"+freq+","+row[7]+","+row[8]+")"
                    self.pointer.execute(insert)
                    counter+=1
        self.pointer.commit()
        self.Disconnect()



    def deleteDB(self):
        os.remove(self.name)











