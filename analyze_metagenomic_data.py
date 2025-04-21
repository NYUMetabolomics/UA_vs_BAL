import os
import sys
import networkx as nx
import math
import statistics
import numpy
numpy.random.seed(19700101)
import random
random.seed(19700101)


heat_threshold = 0.5  # don't show fold changes of magnitude less than ~1.41 
if len(sys.argv) > 1:
    heat_threshold = math.log2(float(sys.argv[1]))

if len(sys.argv) < 3:
    metabolomic_data = None
    for fname in os.listdir("."):
        if fname.endswith("_mx.tsv"):
            metabolomic_data = fname
            break
    if not metabolomic_data:
        print("No file ending with _mx.tsv found! (No metabolomic data found...)")
        sys.exit(-1)
else:
    metabolomic_data = sys.argv[2]

if len(sys.argv) < 4:
    genomic_data = None
    for fname in os.listdir("."):
        if fname.endswith("_gx.tsv"):
            genomic_data = fname
            break
    if not genomic_data:
        print("No file ending with _gx.tsv found! (No genomic data found...)")
        sys.exit(-1)
else:
    genomic_data = sys.argv[3]

network_file = "reaction_network.tsv"
if len(sys.argv) > 4:
    network_file = sys.argv[4]

eliminate_singletons = True  # Don't show disconnected nodes
if len(sys.argv) > 5:
    eliminate_singletons = (sys.argv[5] == "True")

hot_color = "#73FDFF"
if len(sys.argv) > 6:
    hot_color = sys.argv[6]

cold_color = "#FF7E79"
if len(sys.argv) > 7:
    cold_color = sys.argv[7]

wave_number = 0
if len(sys.argv) > 8:
    wave_number = int(sys.argv[8])

transfer_rate = 0.25
if len(sys.argv) > 9:
    transfer_rate = float(sys.argv[9])

node_label = {}
node_fc = {}
def load_nodes(fname):    
    with open(fname, 'r') as file:
        file.readline()  # headers
        for line in file:
            vals = line.strip().split("\t")
            node_id = vals[0]
            fc = float(vals[1])
            label = vals[2]
            if (node_id not in node_fc) or (abs(fc) > abs(node_fc[node_id])):
                node_fc[node_id] = fc
                node_label[node_id] = label


load_nodes(metabolomic_data)
load_nodes(genomic_data)  # We are assuming that genomic data IDs and metabolomic data IDs never overlap (which is safe in the case of InChIK IDs)

G = nx.DiGraph()
with open(network_file, 'r') as myfile:
    for line in myfile:
        vals = line.strip().split("\t")
        if vals[0] == "node":
            node_id = vals[2]
            G.add_node(node_id)
            G.nodes[node_id]["class"] = vals[1]
            G.nodes[node_id]["heat"] = 0
            if node_id in node_fc:
                G.nodes[node_id]["fc"] = node_fc[node_id]
                G.nodes[node_id]["heat"] = abs(node_fc[node_id])
        if vals[0] == "edge":
            source = vals[1]
            dest = vals[2]
            G.add_edge(source, dest)

tot_nodes = len(G.nodes)
tot_edges = len(G.edges)
tot_metabolites = 0
tot_measured_metabolites = 0
tot_reactions = 0   
tot_measured_reactions = 0
metabolite_abs_fcs = []
ortholog_abs_fcs = []

for node in G.nodes:
    if G.nodes[node]["class"] == "metabolite":
        tot_metabolites += 1
        if node in node_fc:
            tot_measured_metabolites += 1
            metabolite_abs_fcs.append(G.nodes[node]["heat"])
    else:
        tot_reactions += 1
        if node in node_fc:
            tot_measured_reactions += 1
            ortholog_abs_fcs.append(G.nodes[node]["heat"])

print("**********************")
print(f"Reference Network Nodes: {tot_nodes}")
print(f"Reference Network Metabolites: {tot_metabolites}")
print(f"Reference Network Orthologs: {tot_reactions}")
print(f"Reference Network Edges: {tot_edges}")
print("**********************")
print(f"Measured Metabolites: {tot_measured_metabolites}")
print(f"Measured Orthologs: {tot_measured_reactions}")
print(f"Median Metabolite Absolute FC: {2**statistics.median(metabolite_abs_fcs):.2f}")
print(f"Median Ortholog Absolute FC: {2**statistics.median(ortholog_abs_fcs):.2f}")

for t in range(wave_number):
    for n in G.nodes:
        ocount = 0
        otot = 0
        o_heat = 0
        for o in G.neighbors(n):
            ocount += 1
            otot += G.nodes[o]["heat"]
        if otot > 0:
            o_heat = otot / ocount
        G.nodes[n]["o_heat"] = o_heat
    for n in G.nodes:
        G.nodes[n]["heat"] += transfer_rate*G.nodes[n]["o_heat"]



to_remove = set()
for n in G.nodes:
    if G.nodes[n]["heat"] <= heat_threshold:
        to_remove.add(n)

for n in to_remove:
    G.remove_node(n)


if eliminate_singletons:
    to_remove = set()
    for n in G.nodes:
        if G.degree[n] == 0:
            to_remove.add(n)
    for n in to_remove:
        G.remove_node(n)

final_nodes = len(G.nodes)
final_edges = len(G.edges)
final_metabolites = 0
final_reactions = 0
for node in G.nodes:
    if G.nodes[node]["class"] == "metabolite":
        final_metabolites += 1
    else:
        final_reactions += 1

print("**********************")
print(f"Final Network Nodes: {len(G.nodes)}")
print(f"Final Network Metabolites: {final_metabolites}")
print(f"Final Network Orthologs: {final_reactions}")
print(f"Final Network Edges: {len(G.edges)}")
print("**********************")
print(f"Percentile of Reference Nodes: {100.0 * float(final_nodes) / float(tot_nodes):.2f}%")
print(f"Percentile of Reference Edges: {100.0 * float(final_edges) / float(tot_edges):.2f}%")
print("**********************")

for node in G.nodes:
    if node not in node_label:
        node_label[node] = node

from matplotlib import colors
from matplotlib.colors import LinearSegmentedColormap

color_list = [colors.to_rgb(cold_color) ,(1,1,1), colors.to_rgb(hot_color)]
cmap = LinearSegmentedColormap.from_list('custom', color_list, N=256)


def fc(fc_in):
    fc_out = (fc_in + 3) / 6.0
    return colors.to_hex(cmap(fc_out))


viz = open(f"heatwave.html", 'w')

print("""
    <html>
    <head>
    <title>Multi-Omic HeatWave Visualization</title>
    <meta name="viewport" content="width=device-width, user-scalable=no, initial-scale=1, maximum-scale=1">
    <script src="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.23.0/cytoscape.min.js"></script>
    <script src="https://cdn.jsdelivr.net/gh/cytoscape/cytoscape.js-pdf-export@main/dist/cytoscape-pdf-export.js"></script>    
    <style>
            body {
                font-family: helvetica;
                font-size: 14px;
            }

            #cy {
                width: 100%;
                height: 100%;
                position: absolute;
                left: 0;
                top: 0;
                z-index: 999;
            }

            h1 {
                opacity: 0.5;
                font-size: 1em;
            }
        </style>

        <script>
            var checkpoint = null;
            var checkpoint_style = [{
                    selector: 'node',
                    labelValign: 'middle',
                    style: {
                        'text-valign': function(ele) {if (ele.data('shape') === "ellipse") { return "top" } else { return "top" } },
                        'width': function(ele) {if (ele.data('shape') === "ellipse") { return 30 } else { return 30 } },
                        'height': function(ele) {if (ele.data('shape') === "ellipse") { return 30 } else { return 30 } },
                        'font-weight': function(ele) { if (ele.data('shape') === "ellipse") { return 'bold' } else { return "normal" } },
                        'font-size': function(ele) { if (ele.data('shape') === "ellipse") { return '16px' } else { return "12px" } },
                        'color': function(ele) { if (ele.data('shape') === "ellipse") { return 'darkblue' } else { return "black" } },
                        'content': 'data(label)',
                        'border-width': 1,
                        //'border-color': 'black',
                        'background-color': function(ele) { return ele.data('bg') },
                        'shape': function(ele) { return ele.data('shape') }
                    }
                },


                {
                    selector: 'edge',
                    style: {
                        'curve-style': 'bezier',
                        'source-arrow-shape': 'triangle',
                        'target-arrow-shape': 'triangle'
                    },
                    css: {
                        'line-color': '#cccbcb' //'#f92411'
                    }
                },

                {
                    selector: ':selected',
                    style: {
                        'background-color': 'purple',
                        'line-color': 'purple',
                        'source-arrow-color': 'purple',
                        'target-arrow-color': 'purple'
                    }
                }
            ];        
            document.addEventListener('DOMContentLoaded', function(){

                var cy = window.cy = cytoscape({
                    container: document.getElementById('cy'),

                    autounselectify: false,
                    
                    boxSelectionEnabled: true,

                    layout: {
                        name: 'cose',
                        nodeOverlap: 1000,
                        animate: false,
                        //idealEdgeLength: 50,
                        nodeDimensionsIncludeLabels: false,
                    },

                    style: checkpoint_style,

                    elements: {
                      nodes: [
""", file=viz)

for n in G.nodes:
    if "fc" in G.nodes[n]:
        color = fc(G.nodes[n]["fc"])
    else:
        color = "#D3D3D3"
    if G.nodes[n]["class"] == "reaction":
        shape = "rectangle"
        link = "https://www.genome.jp/entry/"
    else:
        shape = "ellipse"
        link = "https://www.ebi.ac.uk/chebi/advancedSearchFT.do?searchString="
    # NOTE: It is essential that nodel_labels _not_ contain double-quotes (")
    print(f"""
        {{
          data: {{
            id: "{n}",
            label: "{node_label[n]}",
            bg: "{color}",
            url: "{link}{n}",
            shape: "{shape}"
          }}
        }},
    """, file=viz)
print("""
                      ],
                      edges: [
""", file=viz)

for edge in G.edges:
    print(f"""
    {{
      data: {{
        source: "{edge[0]}",
        target: "{edge[1]}"
      }}
    }}, 
    """, file=viz)

print("""
                      ]
                    }
                });
                cy.on('tap', 'node', function(){
                    try { // your browser may block popups
                        window.open( this.data('url') );
                    } catch(e){ // fall back on url change
                        window.location.href = this.data('url');
                    }
                }); 

            });
        </script>
    </head>

<body>
    <div id="cy"></div>
    <script type="text/javascript">
    document.addEventListener(
        "keydown",
        (event) => {
            const keyName = event.key;

            if (keyName === "f") {
                the_label = prompt("Node Label:");
                var eles = cy.elements('node[label="' + the_label + '"]');
                eles.select();
                cy.center(eles);
                event.preventDefault();
                return;
            }

            if (keyName === "e") {
                cy.elements(":selected").neighborhood().select()
                event.preventDefault();
                return;
            }

            if (keyName === "d") {
                if (event.altKey) {
                    cy.remove(cy.elements(":unselected"));
                    cy.elements(":selected").unselect();
                } else {
                    cy.remove(cy.elements(":selected"));
                }
                event.preventDefault();
                return;
            }

            if (keyName === "Delete") {
                if (event.altKey) {
                    cy.remove(cy.elements(":unselected"));
                    cy.elements(":selected").unselect();
                } else {
                    cy.remove(cy.elements(":selected"));
                }
                event.preventDefault();
                return;
            }


            if (keyName === "D") {
                cy.remove(cy.elements(":unselected"));
                cy.elements(":selected").unselect();
                event.preventDefault();
                return;
            }

            if (keyName === "l") {
                cy.layout({
                    name: 'cose',
                    nodeOverlap: 1000,
                    animate: false,
                    // idealEdgeLength: 64,
                    nodeDimensionsIncludeLabels: true,
                }).run()
                event.preventDefault();
                return;
            }


            if (keyName === "c") {
                if (event.altKey) {
                    if (checkpoint) {
                        cy.json(checkpoint);
                    }
                } else {
                    checkpoint = cy.json()
                    checkpoint.style = checkpoint_style;
                }
                event.preventDefault();
                return;
            }

            if (keyName === "C") {
                if (checkpoint) {
                    cy.json(checkpoint);
                }
                event.preventDefault();
                return;
            }

            if (keyName === "p") {
                cy.pdf({ bg: '#FFF', save: true, full: true })
            }


            if (keyName === "s") {

                const link = document.createElement("a");
                const content = JSON.stringify(cy.json());
                const file = new Blob([content], { type: 'application/json' });
                link.href = URL.createObjectURL(file);
                link.download = "network.json";
                link.click();
                URL.revokeObjectURL(link.href);
            }

            if (keyName === "n") {
                alert('This network has ' + cy.nodes().length + ' nodes and ' + cy.edges().length + ' edges.')
                event.preventDefault();
                return;
            }

            if (keyName === "h") {
                alert('f = find\\ne = expand selection\\nd = delete selection\\ndelete-key = delete selection\\nalt-d = delete non-selected\\nshift-d = delete non-selected\\nalt-delete-key = delete non-selected\\nl = layout\\nc = checkpoint\\nalt-c = reset to checkpoint\\nshift-c = reset to checkpoint\\np = print to pdf\\ns = save to json_file\\ni = import json_file\\nn = show network stats\\nh = help')
                event.preventDefault();
                return;
            }

            if (keyName === "i") {
                var input = document.createElement('input');
                input.type = 'file';

                input.onchange = e => {
                    var file = e.target.files[0];
                    var reader = new FileReader();
                    reader.readAsText(file, 'UTF-8');

                    reader.onload = readerEvent => {
                        var content = readerEvent.target.result;

                        checkpoint = JSON.parse(content);
                        checkpoint.style = checkpoint_style;
                        cy.json(checkpoint);
                    }
                }

                input.click();
            }
        },
        false,
    );
    </script>
</body>

</html>""", file=viz)

    
viz.close()

######################################################


active = open(f"heatwave.tsv", 'w')
print("Class\tNode_ID\tNode_Label\tLog2FC\tHeat", file=active)
for n in G.nodes:
    heat = ""
    log2fc = ""
    if "fc" in G.nodes[n]:
        log2fc = f'{G.nodes[n]["fc"]}'
    heat = f'{G.nodes[n]["heat"]}'
    line = f"{G.nodes[n]['class']}\t{n}\t{node_label[n]}\t{log2fc}\t{heat}"
    print(line, file=active)
active.close()
