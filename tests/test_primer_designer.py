from FluGibson.primer_designer import PrimerDesigner
import os
import networkx as nx

# Get the directory of this test file.
package_directory = os.path.dirname(os.path.abspath(__file__))

os.chdir(os.path.join(package_directory, 'test_data'))

p = PrimerDesigner()
p.read_sequences('victoria_np.fasta')
p.construct_graph()
p.design_assembly_primers()
p.design_junction_sequencing_primers()
p.design_fragment_sequencing_primers()
p.compute_pcr_protocol()

print(p.pcr_protocol)
for n, d in p.graph.nodes(data=True):
    print(n)
    print(d['fragment_sequencing_primers'])
# print(p.graph.nodes(data=True))


def test_construct_graph():
    assert len(p.graph) == 2
    assert len(p.graph.edges()) == 2
    assert isinstance(p.graph, nx.DiGraph)


def test_compute_pcr_protocol():
    pass


def test_design_fragment_sequencing_primers():
    for n, d in p.graph.nodes(data=True):
        frag_seq_primers = d['fragment_sequencing_primers']

        # Check that the forward primers are designed correctly.
        for i, primer in enumerate(frag_seq_primers['fw']):
            if i >= 1:
                assert str(primer) in str(n.seq)
            else:
                assert str(primer) in str(p.graph.predecessors(n)[0].seq)

        # Check that the reverse primers are designed correctly.
        for i, primer in enumerate(frag_seq_primers['re']):
            if i >= 1:
                assert str(primer) in str(n.seq.reverse_complement())

            else:
                assert str(primer) in str(
                    p.graph.successors(n)[0].seq.reverse_complement())


def test_get_fragment_sequencing_primers():
    primers_df = p.get_fragment_sequencing_primers("Vic_NP")
    assert len(primers_df) == 8
