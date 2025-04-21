import os
import sys
import time
import requests
import json
import networkx as nx
import random
import numpy


VERBOSE = False

allow_unmeasurable_reactions = False

random.seed(19700101)
numpy.random.seed(19700101)

if not os.path.isdir('data_cache'):
    os.mkdir('data_cache')


def cached_resource(prefix, url, error_message, suffix, cache_dir="data_cache"):
    cache_fname = None
    for fname in os.listdir(cache_dir):
        if fname.startswith(prefix):
            cache_fname = fname
            break
    if not cache_fname:
        r = requests.get(url)
        if r.status_code != 200:
            print(error_message)
            print(r.text)
            sys.exit(-1)
        else:
            cache_fname = f"{prefix}_{round(time.time())}.{suffix}"
            with open(os.path.join(cache_dir, cache_fname), 'w') as out:
                print(r.text.strip(), file= out)
    return os.path.join(cache_dir, cache_fname)


#
# Get PubChem_Substances
#

pubchem_substances = cached_resource(
    "PubChem_Substances", 
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sourceall/KEGG/xrefs/RegistryID/json",
    "Failed to get KEGG-related substances from PubChem...",
    "json")

with open(pubchem_substances, 'r') as myfile:
    data=myfile.read()
obj = json.loads(data)

sid_from_kegg = {}
kegg_to_sid = {}
for entry in obj['InformationList']['Information']:
    assert entry["SID"] not in sid_from_kegg
    sid_from_kegg[entry["SID"]] = set()
    for kegg_id in entry["RegistryID"]:
        if kegg_id.startswith("D"):
            continue
        sid_from_kegg[entry["SID"]].add(kegg_id)
        assert kegg_id not in kegg_to_sid
        kegg_to_sid[kegg_id] = entry["SID"]

for sid in sid_from_kegg:
    assert len(sid_from_kegg[sid]) < 2

#
# Get PubChem_InChIKeys (not using cached_resource() only because of the two-phase, listkey-based API)
#

pubchem_inchikeys = None
for fname in os.listdir("data_cache"):
    if fname.startswith("PubChem_InChIKeys_"):
        pubchem_inchikeys = os.path.join("data_cache", fname)
        break
if not pubchem_inchikeys:
    r = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sourceall/KEGG/cids/json?list_return=listkey")
    if r.status_code != 200:
        print("Failed to generate listkey or KEGG-related CIDS from PubChem...")
        print(r.text)
        sys.exit(-1)
    else:
        listkey = json.loads(r.text)["IdentifierList"]["ListKey"]
        r = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey/{listkey}/property/inchikey/json")
        if r.status_code != 200:
            print("Failed to retrieve data associated with listkey:", listkey)
            print(r.text)
            sys.exit(-1)
        else:
            pubchem_inchikeys = os.path.join("data_cache", f"PubChem_InChIKeys_{round(time.time())}.json")
            with open(pubchem_inchikeys, 'w') as out:
                print(r.text.strip(), file= out)

cid_to_inchikey = {}
with open(pubchem_inchikeys, 'r') as myfile:
    data=myfile.read()
obj = json.loads(data)

for entry in obj['PropertyTable']['Properties']:
    assert "CID" in entry.keys()
    cid = entry["CID"]
    assert cid not in cid_to_inchikey
    cid_to_inchikey[cid] = entry["InChIKey"][:14]  # only the core graph...

#
# Get PubChem_SID_CID
#

pubchem_substance_to_compound = cached_resource(
    "PubChem_SID_CID", 
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sourceall/KEGG/cids/json",
    "Failed to get KEGG-related substance_compound information from PubChem...",
    "json")

cid_from_sid_from_kegg = {}
sid_to_cid = {}
with open(pubchem_substance_to_compound, 'r') as myfile:
    data=myfile.read()
obj = json.loads(data)

for entry in obj['InformationList']['Information']:
    if "CID" in entry.keys():
        assert len(entry["CID"]) == 1
        cid = entry["CID"][0]
        sid = entry["SID"]
        sid_to_cid[sid] = cid
        if cid not in cid_from_sid_from_kegg:
            cid_from_sid_from_kegg[cid] = set()
        cid_from_sid_from_kegg[cid].add(sid)

#
# Create KEGG_to_InChIK_table
#

kegg_to_inchik = {}
for k in kegg_to_sid:
    s = kegg_to_sid[k]
    if s in sid_to_cid:
        c = sid_to_cid[s]
        inchik = cid_to_inchikey[c]
        kegg_to_inchik[k] = inchik

#
# Get KEGG compound names from KEGG
#

kegg_compounds = cached_resource(
    "KEGG_compounds", 
    "https://rest.kegg.jp/list/compound",
    "Failed to get KEGG compound names...",
    "tsv")

#
# Find best compound per compound name (compound with the most reactions)
#

compounds = {}
compound_choices = {}
with open(kegg_compounds) as f:
    for line in f:
        # "C00018\tPyridoxal phosphate; Pyridoxal 5'-phosphate; Pyridoxal 5-phosphate; PLP"
        (cname, names) = line.strip().split("\t")
        if ";" in names:
            name = names.split(";")[0]
        else:
            name = names
        if name in compounds:
            if name not in compound_choices:  # This means there is more than one KEGG compound for the same name which _can_ happen!!!
                compound_choices[name] = set([compounds[name]])
            compound_choices[name].add(cname)
        compounds[name] = cname

for name in compound_choices:
    best_rcount = -1
    for cname in sorted(compound_choices[name]):
        #
        # Get information about candidate cname for compound name...
        #
        cname_file = cached_resource(
            cname,
            "https://rest.kegg.jp/get/" + cname,
            f"Failed to get info about {cname}...",
            "txt")
        with open(cname_file) as f:                
            in_reactions = False
            rcount = 0
            for line in f:
                if in_reactions:
                    if line[0] != " ":
                        in_reactions = False
                        continue
                    else:
                        rcount += len(line.strip().split())
                else:
                    if line.startswith("REACTION"):
                        rcount += len(line.split()[1:])
            if rcount > best_rcount:
                best_rcount = rcount
                compounds[name] = cname


#
# Get KEGG reaction-to-ortholog list
#

kegg_reaction_to_ortho = cached_resource(
    "KEGG_rn_to_ko", 
    "https://rest.kegg.jp/link/ko/rn",
    "Failed to get KEGG reaction-to-ortholog list...",
    "tsv")

reaction_to_ortho = {}
with open(kegg_reaction_to_ortho) as f:
    for line in f:
        vals = line.strip().split("\t")
        reaction = vals[0][3:]
        ortholog = vals[1][3:]
        if reaction not in reaction_to_ortho:
            reaction_to_ortho[reaction] = set()
        reaction_to_ortho[reaction].add(ortholog)

reactions = list(reaction_to_ortho.keys())

#
# Get KEGG glycan names from KEGG  <-- sometimes they are referenced by name in reactions rather than as GXXXXX
#

kegg_glycans = cached_resource(
    "KEGG_glycans", 
    "https://rest.kegg.jp/list/gl",
    "Failed to get KEGG glycan names...",
    "tsv")

glycans = {}
with open(kegg_glycans) as f:
    for line in f:
        # "G00001   N-Acetyl-D-glucosaminyldiphosphodolichol; (GlcNAc)1 (PP-Dol)"
        (glname, names) = line.strip().split("\t")
        if ";" in names:
            name = names.split(";")[0]
        else:
            name = names
        if name in glycans:
            if VERBOSE:
                print("Note: recurring glycan name (linear representation of tree will do this)...", name, glycans[name], glname)
            continue
        if name in compounds:
            if VERBOSE:
                print("Warning: name shared by glycans and compounds!!!", name, compounds[name], glname)
            continue  # we will assume that if a reaction references this name we can use the compound rather than the glycan
        glycans[name] = glname

#
# Get KEGG reactions
#

kegg_reactions = cached_resource(
    "KEGG_reactions", 
    "https://rest.kegg.jp/list/rn",
    "Failed to get KEGG reactions...",
    "tsv")


def parse_reaction_side(reaction_side):
    inchikeys = []
    lparts = reaction_side.split(" + ")
    lsplinters = []
    for lpart in lparts:
        if lpart.endswith("(side 1)") or lpart.endswith("(side 2)"):
            lpart = lpart[:lpart.rfind("(")]  # this must be done exceptionally becase of the space character inside the parens...
        lsplinters = lpart.split(" ")
        if lsplinters[0].isdigit() or \
           lsplinters[0] == "n" or \
           lsplinters[0] == "n-1" or \
           lsplinters[0] == "(n+1)" or \
           lsplinters[0] == "(n-1)" or \
           lsplinters[0] == "(n-2)" or \
           lsplinters[0] == "2n" or \
           lsplinters[0] == "4n":
            lsplinters = lsplinters[1:]

        no_quant_suffix = []

        for lsplinter in lsplinters:
            if lsplinter.endswith("(n)") or \
               lsplinter.endswith("(x)") or \
               lsplinter.endswith("(m-1)") or \
               lsplinter.endswith("(m+1)") or \
               lsplinter.endswith("(n-1)") or \
               lsplinter.endswith("(n-x)") or \
               lsplinter.endswith("(n+1)") or \
               lsplinter.endswith("(n+2)") or \
               lsplinter.endswith("(m)") or \
               lsplinter.endswith("(n+m)") or \
               lsplinter.endswith("(m+n)") :
                lsplinter = lsplinter[:lsplinter.rfind("(")]
            assert lsplinter  # if lsplinter == "" this means we need to add a case to the lsplinter[0] if statement above...
            no_quant_suffix.append(lsplinter)

        lsplinters = no_quant_suffix

        lpart = " ".join(lsplinters)

        if lpart.startswith("G") and lpart[1:].isdigit():
            continue  # ignore GXXXXX substance in reaction

        if lpart in glycans:
            continue

        if lpart not in compounds:
            print("Unrecognized compound", lpart)
            print(f"Found in reaction {rname}: {reaction}")
            sys.exit(-1)

        cname = compounds[lpart]
        if cname not in kegg_to_inchik:
            if VERBOSE:
                print("Cannot assign InChiK to KEGG compound:", cname)
                print("Which is the KEGG CXXXXX ID for:", lpart)
                print(f"Found in reaction {rname}: {reaction}")
        else:
            inchikeys.append(kegg_to_inchik[cname])
    return inchikeys


inchik_left = {}
inchik_right = {}

with open(kegg_reactions) as f:
    for line in f:
        # "R00001\tpolyphosphate polyphosphohydrolase; Polyphosphate + n H2O <=> (n+1) Oligophosphate"
        (rname, details) = line.strip().split("\t")

        if rname not in reactions:
            if allow_unmeasurable_reactions:
                reactions.add(rname)
            else:
                continue

        if ";" in details:
            last_semicolon = details.rfind(";")
            name = details[:last_semicolon]
            reaction = details[(last_semicolon + 2):]  # The +1 is because of a space character after the semicolon
        else:
            reaction = details
        (left, right) = reaction.split(" <=> ")

        left_inchiks = parse_reaction_side(left)
        right_inchiks = parse_reaction_side(right)

        for inchik in left_inchiks:
            if rname not in inchik_left:
                inchik_left[rname] = set()
            inchik_left[rname].add(inchik)

        for inchik in right_inchiks:
            if rname not in inchik_right:
                inchik_right[rname] = set()
            inchik_right[rname].add(inchik)

sorted_reactions = sorted(reactions)

G = nx.DiGraph()
for reaction in sorted_reactions:
    if (reaction not in inchik_right) and (reaction not in inchik_left):  # we ignore "standalone" reactions...
        continue
    representatives = set()
    for gene in reaction_to_ortho[reaction]:
        if gene not in G.nodes:
            G.add_node(gene)
            G.nodes[gene]["class"] = "reaction"
        representatives.add(gene)
    # else:
    #     G.add_node(reaction)
    #     G.nodes[reaction]["class"] = "reaction"
    #     representatives.add(reaction)  # We are making the assumption no gene was ever called R#####...
    if reaction in inchik_left:
        for inchik in inchik_left[reaction]:
            if inchik not in G.nodes:
                G.add_node(inchik)
                G.nodes[inchik]["class"] = "metabolite"
            for rep in representatives:
                G.add_edge(inchik, rep)
    if reaction in inchik_right:
        for inchik in inchik_right[reaction]:
            if inchik not in G.nodes:
                G.add_node(inchik)
                G.nodes[inchik]["class"] = "metabolite"
            for rep in representatives:
                G.add_edge(rep, inchik)

U = G.to_undirected()

with open(f"reaction_network.tsv", 'w') as out:
    # print("EntryType\tNode_ID_or_From_ID\tNode_Label_or_To_ID", file=out)
    for n in U.nodes:
        line = f"node\t{U.nodes[n]['class']}\t{n}"
        print(line, file=out)
    for (a, b) in U.edges:
        print(f"edge\t{a}\t{b}", file=out)
