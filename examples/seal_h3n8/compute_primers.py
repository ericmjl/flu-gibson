import os
import pandas as pd
import sys

from FluGibson.primer_designer import PrimerDesigner
from FluGibson.plasmids import plasmids

segmaps = dict()
segmaps['seg1'] = 'PB2'
segmaps['seg2'] = 'PB1'
segmaps['seg3'] = 'PA'
segmaps['seg4'] = 'HA'
segmaps['seg5'] = 'NP'
segmaps['seg6'] = 'NA'
segmaps['seg7'] = 'M'
segmaps['seg8'] = 'NS'

if __name__ == '__main__':

    master_cloning_primers = list()

    for f in os.listdir(os.getcwd()):
        if '.fasta' in f:
            print(f)
            for k, v in segmaps.items():
                if k in f:
                    segment = v
            backbone = plasmids['pDZ-{0}'.format(segment)]
            p = PrimerDesigner()
            p.read_sequences(f)
            p.add_sequences([backbone])
            p.construct_graph()
            protocol = p.pcr_protocol()

            f = f.replace('.fasta', '')
            pd.DataFrame(protocol).to_csv('{0}_pcr_protocol.csv'.format(f))
            for n, d in p.nodes(data=True):
                # print(n, d)
                # for direction in ['fw', 're']:
                #     primer_data = dict()
                #     primer_data['{0}_junction_primer'.format(direction)] = \
                #         d['{0}_sequence'.format(direction)].seq
                #     primer_data['{0}_sequencing_primer'.format(direction)] = \
                #         d['{0}_sequencing_primer'.format(direction)].seq

                #     primer_data['part'] = n.id
                #     primer_data['{0}_junction_primer_name'.format(direction)] = '{0}_{1}_junction'.format(direction, n.id)
                #     primer_data['{0}_seuqencing_primer_name'.format(direction)] = '{0}_{1}_sequencing'.format(direction, n.id)

                #     master_primer_list.append(primer_data)

                """Grab out the cloning primers."""
                primer_data = dict()
                for drxn in ['fw', 're']:
                    primer_data['part'] = n.id
                    primer_data['{0}_cloning_primer'.format(drxn)] = str(d['{0}_cloning_primer'.format(drxn)].seq)
                master_cloning_primers.append(primer_data)

    pd.DataFrame(master_cloning_primers).to_csv('all_primers.csv')
