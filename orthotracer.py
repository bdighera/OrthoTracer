import ast, os
from Bio import SeqIO
from time import sleep

from utils.visualizations import PhyloTreeConstruction
from utils.argparse import cli
from utils.collector import SequenceCollector, GenomicContext
from dB.sqlite import SQliteRecordInput, VerifyDatabase, RecordRetrival


'''
This is the main file for the second module of the pipeline - retrieve records from records.db and build figure.
Pipeline has been figured to call Main() from juypter-notebook with JupyterNotebook=True
'''

class orthotracer():

    def __init__(self):
        pass

    def visualizeRecords(self, accessionList = None, fileHandle = '', prune='', JupyterNotebook= False, printLeafAccession= False,
         printLeafDesc= False,
         image_scaling = 1,
         stretch= 0,
         inputFile = os.path.join('samples', 'Human_HSP70_paralogs.txt')):

        '''
        :param accessionList: - fasta accession list file path as input (Experimental - not sure if working)
        :param fileHandle: file that should be used for input - not fasta fomatted, just accession numbers line after line
        :param prune: list of sequences, will remove from figures
        :param JupyterNotebook: Toggle for jupyter-notebook compatibility
        :param printLeafAccession: Toggle the printing of the leaf accessions
        :param printLeafDesc: Toggle the printing of the leaf descriptions
        :param image_scaling: used in domain figure - increase if the motifs are far away from the phylogenetic tree
        :param stretch: used in genomic context figure - increase if the motifs are far away from the phylogenetic tree
        :return: Figures
        '''

        JupyterNotebookInlineFig = JupyterNotebook

        V = RecordRetrival()

        if JupyterNotebookInlineFig == True:
            if accessionList != None:

                V.retrieveFastaRecords(accessionList)

            elif accessionList == None and fileHandle == None:
                raise ValueError('Need to pass list of accession numbers to Main.py')

            elif fileHandle != '':
                V.retrieveFastaRecords(fileHandle)
            pass

        elif JupyterNotebookInlineFig == False:
            V.retrieveFastaRecords(inputFile)

        records = V.pullDBrecords(dbfile='dB/Records.db')

        proteinAccession = [record[0][1] for record in records]
        proteinSeq = [record[0][2] for record in records]
        proteinDescription = [record[0][3] for record in records]
        geneID = [record[0][11] for record in records]
        genomicContext = [ast.literal_eval(record[0][12]) for record in records]
        parentDomains = [ast.literal_eval(record[0][13]) for record in records]
        introns = [record[0][14] for record in records]
        exonLength = [record[0][15] for record in records]
        taxonomy = [record[0][16] for record in records]
        commonNames = [record[0][17] for record in records]

        Phylo = PhyloTreeConstruction(proteinAccession,
                                      proteinSeq,
                                      proteinDescription,
                                      genomicContext,
                                      parentDomains,
                                      introns,
                                      exonLength,
                                      JupyterNotebookInlineFig,
                                      prune,
                                      printLeafAccession,
                                      printLeafDesc,
                                      commonNames,
                                      geneID,
                                      image_scaling,
                                      stretch)

        return Phylo


    def scrapeRecords(self, path='input.txt', db='dB/Records.db', linux=True, GCsize=50000, apikey='4e3f380c489dcaacecf12c2c3483ebe24909'):

        # Args = ArguementParser()
        #
        # proteinPath, cwd, rec, outfile, genomicPath, CDSPath, tree = Args.parse()

        # sets the timer value between ping requests to the NCBI server
        timer = 6

        proteinSeqRecords = SeqIO.parse(path, 'fasta')

        # Iterates through each parsed record in the input file
        for proteinRecord in proteinSeqRecords:

            try:

                # Uses the record accession number to check if it is already in the SQL database, if not: continue.
                SQL = VerifyDatabase(proteinRecord.id, db=db)

                # If accession number is not located in the DB then it will run through the pipeline
                if SQL.checkRecords() != True:

                    # Initialize Sequence Collector class with input of the working protein accession record
                    Seq = SequenceCollector(proteinRecord, linux=linux, apikey=apikey)

                    # collects the protein ID for the current working protein
                    # protein ID is needed in order to establish relationship between NCBI databases ie. Gene, Protein, etc.
                    # runs timer after the method is executed in order to not ping NCBI server too fast - or they will SHUT YOU DOWN
                    id = Seq.collectProteinIDs()
                    sleep(timer)

                    # Runs the method to collect the CDS and taxonomy
                    # TO understand this method better (along with other methods of class Sequence Collector) visit SequenceCOllector.py in CWD
                    CDS, Taxonomy = Seq.collectCDS()
                    sleep(timer)

                    # Collects the entire sequence of the gene and the gene ID
                    Genomic, GeneID = Seq.collectGenomicSeq()
                    sleep(timer)

                    # Collects the intron phases and exon lengths - this method requires executible files from ./exec folder
                    # Internal function -does not ping NCBI server
                    IntronPhases, ExonLengths = Seq.collectIntrons(CDS, Genomic)

                    # Collects the domains for the protein
                    Domains = Seq.collectProteinDomains()
                    sleep(timer)

                    # Initialize Genomic Context class
                    # Collects neighboring genes to parent protein, coding direction, and domains of those genes
                    GC = GenomicContext(Genomic, kbps=GCsize, apikey=apikey)

                    # Collects record for  +/- 50k basepairs up/downstream of parent gene
                    gcRecord = GC.fetchRecord()
                    sleep(timer)

                    # parses raw genomic context data collected from NCBI
                    parsedGCrecord = GC.parseRecord(gcRecord)

                    commonName = Seq.collectCommonName(id)

                    # Initialize the SQLite class with all of the data for the working protein
                    SQL = SQliteRecordInput(Seq.proteinRecord,
                                            Seq.proteinID,
                                            CDS,
                                            Genomic,
                                            GeneID,
                                            parsedGCrecord,
                                            Domains,
                                            IntronPhases,
                                            ExonLengths,
                                            Taxonomy,
                                            commonName,
                                            db)

                    # Uploads records to the DB
                    # REMEMBER TO CHANGE THE NAME OF THE TABLE TO REFLECT WHERE YOU WANT THE RECORDS TOGO!
                    SQL.uploadRecords()

                # If protein accession number is in the database then it will skip over and go to the next record
                else:
                    print('Record ' + proteinRecord.id + ' is already in database. Proceeding to next record')

            # Errors in some of the sequences do occur on the NCBI server end. For this errors the DNAJC_Errors file is populated with those accession #s
            except IndexError:
                print('Index Error Occurred on: ' + str(proteinRecord.id))
                with open('DNAJC_Errors', 'a') as handle:
                    handle.write(str(proteinRecord.id) + '\n')
                    handle.close()


if __name__ == '__main__':
    cli = cli()
    ot = orthotracer()
    print(cli)

    if cli['c'] != False:
        ot.scrapeRecords(path=cli['f'],
                          db=cli['db'],
                          linux=cli['l'],
                          GCsize=cli['s'],
                          apikey=cli['apikey'])

    elif cli['v'] != False:
        records = ot.visualizeRecords(inputFile=cli['f'])
        if cli['i'] == True:
            records.renderingTreewithIntrons(saveFile=cli['saveFig'], showFig=cli['showFig'])
        if cli['gc'] == True:
            records.renderingTreewithGenomicContext(saveFile=cli['saveFig'], showFig=cli['showFig'])
        if cli['d'] == True:
            records.renderingTreewithDomains(saveFile=cli['saveFig'], showFig=cli['showFig'])




