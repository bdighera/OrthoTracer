import argparse


def cli():

        # Initializes the parsing function
        parser = argparse.ArgumentParser()
        group = parser.add_mutually_exclusive_group(required=True)


        parser.add_argument('-f', '-file', default=None, type=str, help='input file path', required=True)

        group.add_argument('-c','-collect', help='collects protein sequences and transfers to database', action='store_true')
        parser.add_argument('-dB', default='dB/Records.db', type=str, help='path to database, default ./dB/Records.db')
        parser.add_argument('-l', '-linux', default=True, type=bool, help='Set to False if running on MacOX. Warning: Database will populate N/A intron phase on MacOX')
        parser.add_argument('-s', '-size', default=50000, type=int, help='Size of genomic region to collect in kbps')
        parser.add_argument('-apikey', default='4e3f380c489dcaacecf12c2c3483ebe24909', type=str, help='NCBI API key')




        group.add_argument('-v','-visualize', help='visualizes data, ensure that sequences have first been collected', action='store_true')
        parser.add_argument('-i', '-introns', default=True, type=bool, help='switch to turn on and off intron figure. Default displays figure')
        parser.add_argument('-gc', default=True, type=bool, help='switch to turn on and off intron figure. Default displays figure')
        parser.add_argument('-d', '-domains', default=True, type=bool, help='switch to turn on and off domain figure. Default displays figure')
        parser.add_argument('-showFig', default=False, type=bool, help='switch to turn on and off figure saving. Default does not save figure')
        parser.add_argument('-viewfig', default=True, type=bool, help='switch to turn on and off figure viewing. Default displays figure')




        args = vars(parser.parse_args())

        return args

