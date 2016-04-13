from FluGibson.primer_designer import PrimerDesigner

import os
import pandas as pd

master_primer_list = list()

for f in os.listdir(os.getcwd()):
    if '.fasta' in f:
        p = PrimerDesigner()
        p.read_sequences(f)
        p.construct_graph()
        p.compute_pcr_protocol()
        p.save_pcr_protocol()

        graph = p.graph

        for n, d in graph.nodes(data=True):
            print(n, d)
            for direction in ['fw', 're']:
                primer_data = dict()
                primer_data['{0}_junction_primer'.format(direction)] = d['{0}_sequence'.format(direction)].seq
                primer_data['{0}_sequencing_primer'.format(direction)] = d['{0}_sequencing_primer'.format(direction)].seq

                primer_data['part'] = n.id
                primer_data['{0}_junction_primer_name'.format(direction)] = '{0}_{1}_junction'.format(direction, n.id)
                primer_data['{0}_seuqencing_primer_name'.format(direction)] = '{0}_{1}_sequencing'.format(direction, n.id)

                master_primer_list.append(primer_data)

pd.DataFrame(master_primer_list).to_csv('all_primers.csv')
