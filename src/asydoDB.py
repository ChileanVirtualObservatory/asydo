import sqlite3 as lite
import sys
from astropy.io.votable.tree import Field as pField
import sys
import urllib
import urllib2
import requests
import csv

SqlEquivalent = {
    "char":     "TEXT",
    "double":   "DOUBLE",
    "int":      "INT",
    "boolean":  "BOOLEAN"
}


class SplataDBManager:

    name = ""
    connected = False
    fields = {}
    pointer = None
    slap_serv = 'https://find.nrao.edu/splata-slap/slap'
    alma_band_freq = {'3': [88, 116], '4': [125, 163], '6': [211, 275], '7': [275, 373], '8': [385, 500], '9': [602, 720]}

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

    def VOGetLines(log, prefix):
        c = 299792458.0
        for band in alma_band_freq:
            log.write('band=' + band + '\n')
            w_init = c / (alma_band_freq[band][0] * 1000000000.0)
            w_end = c / (alma_band_freq[band][1] * 1000000000.0)
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

    def VOloadFields(self, fields):
        # recieves the fields from the VOtable and adds them to the class
        for f in fields:
            if isinstance(f, pField):
                name = f.name
                description = f.description
                type = SqlEquivalent[f.datatype]
                self.fields[name] = (description, type)

    def createDBFromVO(self):
        command = "CREATE TABLE Catalogo (ID INT PRIMARY KEY NOT NULL,"
        metadata = "CREATE TABLE Metadata (ID INT PRIMARY KEY NOT NULL,"
        metadata += "Column TEXT NOT NULL, Description TEXT NOT NULL)"
        insertmetadata = []
        contador = 0
        for i in self.fields:
            name = i

            descripcion, tipodato = self.fields[i]
            if contador != 0:
                command = command + ", "

            command = command + " " + name.replace(" ", "_") + " " + tipodato
            command2 = "INSERT INTO Metadata VALUES(" + str(
                contador) + ", '" + name + "', '" + descripcion + "')"
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

    def createDBFromCSV(self, filename,output):
        self.Connect(output)
        create = "CREATE TABLE Lines(ID INT PRIMARY KEY NOT NULL,SPECIES TEXT,CHEM_NAME TEXT,FREQ REAL,INTENSITY REAL,EL REAL)"
        drop = "DROP TABLE Lines"
        with open(filename, 'rb') as csvfile:
            sreader = csv.reader(csvfile, delimiter=':', quotechar='|')
            counter=0
            for row in sreader:
                if counter == 0:
                    self.pointer.execute(drop)
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


mng=SplataDBManager()
mng.createDBFromCSV('bin/splatalogue.csv','lines.db')


