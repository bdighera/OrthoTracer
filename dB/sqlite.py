import sqlite3, re
from uuid import uuid4

'''
These classes control of the SQLite DB. Remember that there must be manual editing for the SQL table names
This is needed to make sure that the accessions are checked against/added to the correct table in the DB.
If you require a new table in the db then run this script and update the method makeNewTable of class makeNewTable to reflect the necessary table name
Remember to update that in the other two classes as well
'''

#Takes all of the records from the various collection procedures and formats them, and inputs them into the database
class SQliteRecordInput():

    def __init__(self, Protein, ProteinID, CDS, Genomic, GeneID, GC, Domains, IntronPhase, ExonLength, Taxonomy, CommonName, db):
        self.conn = sqlite3.connect(db)
        self.ProteinRecord = Protein
        self.ProteinSeq = str(Protein.seq)
        self.ProteinAccession = Protein.id
        self.ProteinDescription =Protein.description
        self.ProteinID = int(ProteinID)
        self.CDS = CDS
        self.CDSSeq = str(CDS.seq)
        self.CDSAccession = CDS.id
        self.CDSDescription = CDS.description
        self.Genomic = Genomic
        self.GenomicSeq = str(Genomic.seq)
        self.GenomicAccession = str(Genomic.id)
        self.GenomicDescription = Genomic.description
        self.GeneID = int(GeneID)
        self.GC = str(GC)
        self.Domains = str(Domains)
        self.Introns = str(IntronPhase)
        self.ExonLength = str(ExonLength)
        self.uuid = str(uuid4())
        self.taxonomy = str(Taxonomy)
        self.CommonName = CommonName

    def uploadRecords(self):

        C = self.conn.cursor()

        try:
            C.execute('''INSERT INTO Sequences VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''',(self.uuid, self.ProteinAccession, self.ProteinSeq, self.ProteinDescription, self.ProteinID,
                                                                                  self.CDSAccession, self.CDSSeq, self.CDSDescription, self.GenomicAccession,
                                                                                  self.GenomicSeq, self.GenomicDescription,
                                                                                  self.GeneID,
                                                                                  self.GC,
                                                                                    self.Domains,
                                                                                  self.Introns,
                                                                                  self.ExonLength,
                                                                                    self.taxonomy,
                                                                                          self.CommonName))
        except sqlite3.IntegrityError:
            pass
        self.conn.commit()
        self.conn.close()

#Checks to make sure that the database does not already contain a record before running the NCBI mining procedures
class VerifyDatabase():
    def __init__(self, proteinAccession, db):
        self.connect = sqlite3.connect(db)
        self.proteinAccession = proteinAccession

    #Checks to make sure that the record isn't already in the table before it runs it. Checks against the accession number of the protein fasta
    def checkRecords(self):
        C = self.connect.cursor()

        C.execute('SELECT ProteinAccession FROM Sequences WHERE ProteinAccession= (?)''',(self.proteinAccession,))

        data = C.fetchall()

        if not data:
            return False
        else:
            return True

#Pulls records from the database based on the input sequence list
class PullSequences():
            def __init__(self, proteinAccession, dbfile):
                self.dbfile = dbfile
                self.connect = sqlite3.connect(dbfile)
                self.proteinAccessionList = proteinAccession

            # Checks to make sure that the record isn't already in the table before it runs it. Checks against the accession number of the protein fasta
            def checkRecords(self):

                dataList = []

                C = self.connect.cursor()

                for accession in self.proteinAccessionList:

                    C.execute('SELECT * FROM Sequences WHERE ProteinAccession= (?)''', (accession,))

                    data = C.fetchall()

                    if data != []:
                        dataList.append(data)

                return dataList


#Checks input file for formatting, pulls sequences in for use by SQLiteChecker
class RecordRetrival():

    def __init__(self):
        pass

    # This method is used to retrieve desired records for visualization from input file
    # input file consists of file in cd with fasta formatting
    def accessionExtractor(self, handle):
        print()
        p = re.compile('[NX]P_\d+.\d')
        accessionList = p.findall(str(handle))
        return accessionList

    def retrieveFileRecords(self, filename):
        with open(filename, 'r') as handle:
            contents = handle.readlines()
            contents = [i.strip('\n') for i in contents]

            handle.close()
        self.inputAccessionList = contents

        return self.inputAccessionList

    def retrieveFastaRecords(self, filename):
        with open(filename, 'r') as handle:
            content = handle.read()
        self.inputAccessionList = self.accessionExtractor(content)

        return self.inputAccessionList

    # This method is used to take the input accession list and use it as primary key for sqlite db
    # pulls the db entry from the accession number of input file
    def pullDBrecords(self, dbfile):
        inputList = self.inputAccessionList

        # SQLiteChecker class takes the input accession list and the name of the dbfile
        # dbfile currently set to records.db - only current records file
        SQL = PullSequences(proteinAccession=inputList, dbfile=dbfile)

        records = SQL.checkRecords()

        return records

    def retrieveRecordsbyList(self, accessionList):
        self.inputAccessionList = accessionList
        return self.inputAccessionList

#Makes a brand new table in the database so that different sets of proteins can be added
class makeNewTable():
    def __init__(self):
        self.connect = sqlite3.connect('Records.db')
    def makeNewTable(self):
        C = self.connect.cursor()

        '''
        Insert the desired table name where it says insert table name
        '''
        C.execute('''CREATE TABLE InsertTableNameHere (UUID TEXT PRIMARY KEY ,ProteinAccession TEXT UNIQUE, ProteinSequence TEXT, ProteinDescription TEXT,
                                  ProteinID INTEGER, CDSAccession TEXT, CDSSeq TEXT, CDSDescription TEXT, GenomicAccession TEXT,
                                GenomicSeq TEXT, GenomicDescription TEXT, GeneID INTEGER, GenomicContext TEXT, ParentDomains TEXT,
                                Introns TEXT, ExonLength TEXT, Taxonomy TEXT, CommonName TEXT)''')

        self.connect.commit()

if __name__ == '__main__':

    SQL = makeNewTable()
    SQL.makeNewTable()