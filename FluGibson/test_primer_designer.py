from primer_designer import PrimerDesigner
import os
import networkx as nx

package_directory = os.path.dirname(os.path.abspath(__file__))  # Get the directory of this test file.

os.chdir(os.path.join(package_directory, 'test_data'))

p = PrimerDesigner('victoria_np.fasta')
p.construct_graph()
p.design_primers()
p.compute_pcr_protocol()
p.save_pcr_protocol()

print(p.pcr_protocol)
print(p.graph)

def test_construct_graph():
	assert len(p.graph) == 2
	assert len(p.graph.edges()) == 2
	assert isinstance(p.graph, nx.DiGraph)

def test_compute_pcr_protocol():
	pass