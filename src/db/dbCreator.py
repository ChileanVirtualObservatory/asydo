import sys
import time
import urllib
import os

from asydoDB import DataBase

default_url = "http://www.csrg.cl/~maray/splatalogue.csv"
default_csv_name = "lines2.csv"
default_db_name = "ASYDO"
csv = False
URI = ""
log = sys.stdout

database = DataBase(default_db_name)

def reporthook(count, block_size, total_size):
    global start_time
    if count == 0:
        start_time = time.time()
        return
    duration = time.time() - start_time
    progress_size = int(count * block_size)
    speed = int(progress_size / (1024 * duration))
    percent = int(count * block_size * 100 / total_size)
    log.write("\r\t...%d%%, %d MB, %d KB/s, %d seconds passed" % (percent, progress_size / (1024 * 1024), speed, duration))
    log.flush()



if (len(sys.argv)==1):

    log.write("No source was defined for the Line Database, so the Default Line Database will be imported: \n\tDownloading Default Line CSV \n")
    urllib.urlretrieve(default_url,default_csv_name , reporthook)
    URI = default_csv_name
    log.write("\n\tDefault Line CSV was Downloaded\n")
    csv = True

elif (len(sys.argv)>2):
    URI = sys.argv[2]

    if (sys.argv[1]=="-C"):
        log.write("Using %s as source for Line Database\n" % os.path.basename(URI))
        csv = True

    elif sys.argv[1]=="-T":
        w_range = [88,720]
        message = "No range specified, using Default range (88 Ghz to 720 Ghz)\n"
        if len(sys.argv) == 5 and sys.argv[3] == "-R":
            w_range = sys.argv[4].split(":")
            message = "The specified range is ("+w_range[0]+" to " + w_range[1]+")\n"

        log.write(message)
        database.deleteDB()
        database.VOGetLines(log,URI,w_range)
        database.loadVoTable("./votables/customVOTable.xml",{3:"FREQ",4:"SPECIES",5:"CHEM_NAME",7:"INTENSITY",11:"EL"})




    #elif (sys.argv[1]=="-s"): TODO Allow use of SQLITE (CDMS)

    #else TODO Include helper -h and -help


if csv:
    log.write("Importing CSV (%s) to SQL Database\n" % os.path.basename(URI))
    database.createDBFromCSV(URI,log)

log.write("Database creation is now complete, you can now use ASYDO. Have a nice day.\n")