from FluGibson.primer_designer import PrimerDesigner

import os
import pandas as pd

master_primer_list = list()

for f in os.listdir(os.getcwd()):
    if '.fasta' in f:
        pdesigner = PrimerDesigner(f)
        pdesigner.construct_graph()
        pdesigner.design_primers()
        pdesigner.compute_pcr_protocol()
        pdesigner.save_pcr_protocol()

        graph = pdesigner.graph

        for n, d in graph.nodes(data=True):
            for prefix in ['fw', 're']:
                primer_data = dict()
                primer_data['primer'] = d[
                    '{prefix}_primer'.format(prefix=prefix)]
                primer_data['sequence'] = d[
                    '{prefix}_sequence'.format(prefix=prefix)].seq

                master_primer_list.append(primer_data)

pd.DataFrame(master_primer_list).to_csv('all_primers.csv')