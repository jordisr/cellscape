'''
Read domains from Uniprot XML into a custom object
'''

import xml.etree.ElementTree as ET
import sys
import argparse
import json

class protein:
    """Data structure to hold topological/domain information"""
    def __init__(self,name):
        self.name = name
        self.domains = []
        self.topology = []
        self.ptm = {}
        self.sequence = ""

    def add_domain(self,name, start, end):
        self.domains.append((name, int(start), int(end)))
    def add_topology(self, name, start, end):
        self.topology.append((name, int(start), int(end)))
    def add_ptm(self, name, start, end):
        self.ptm[name] = (int(start), int(end))
    def process_segments(self):
        if 'chain' in self.ptm:
            (self.chain_start, self.chain_end) = self.ptm['chain']
        else:
            (self.chain_start, self.chain_end) = (1, 99999)

        last = self.chain_start
        self.domain_segments = []

        for domain in self.domains:
            if (domain[1] - last) > 1:
                self.domain_segments.append(('None',last, domain[1]-1))

            self.domain_segments.append(domain)
            last = domain[2]

        if (self.chain_end - last) > 1:
            self.domain_segments.append(('None',last, self.chain_end))

def parse_xml(xmlpath):
    tree = ET.parse(xmlpath)
    root = tree.getroot()
    ns = '{http://uniprot.org/uniprot}'
    sequences = []

    for entry in tree.iter(tag=ns+'entry'):
        accession = entry.find(ns+'accession')
        sequence = protein(accession.text)

        for feature in entry.iter(tag=ns+'feature'):

            # look for transmembrane regions
            if feature.get('type') in ('topological domain','transmembrane region'):
                begin = feature.find(ns+'location').find(ns+'begin').get('position')
                end = feature.find(ns+'location').find(ns+'end').get('position')
                feature_description = feature.get('description').split(';')[0]
                sequence.add_topology(feature_description, begin, end)

            # look for protein domains
            elif feature.get('type') == 'domain':
                begin = feature.find(ns+'location').find(ns+'begin').get('position')
                end = feature.find(ns+'location').find(ns+'end').get('position')
                sequence.add_domain(feature.get('description'),begin,end)

            # look for signal peptide and mature chain
            elif feature.get('type') in ('chain', 'propeptide','signal peptide'):
                begin = feature.find(ns+'location').find(ns+'begin').get('position')
                end = feature.find(ns+'location').find(ns+'end').get('position')
                sequence.add_ptm(feature.get('type'), begin, end)

        sequence.process_segments()
        sequences.append(sequence)

        for seq in entry.iter(tag=ns+'sequence'):
            if seq.text is not None:
                sequence.sequence = seq.text.replace('\n','')

    return(sequences)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Parse UniProt XML file',  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--xml', help='Input XML file', required=True)
    parser.add_argument('--json', action='store_true', default=False, help='Output relevant information in JSON')
    args = parser.parse_args()

    uniprot = parse_xml(args.xml)

    for entry in uniprot:
        if args.json:
            data = {
            'name': entry.name,
            'sequence': entry.sequence,
            'domains': entry.domain_segments,
            'topology': entry.topology
            }
            print(json.dumps(data, indent=2))
        else:
            with open(entry.name+'.domains.csv','w') as f:
                f.write(','.join(['res_start','res_end','description'])+'\n')
                for domain in entry.domain_segments:
                    f.write(','.join(map(str,[domain[1],domain[2],domain[0]]))+'\n')

            with open(entry.name+'.topology.csv','w') as f:
                f.write(','.join(['res_start','res_end','description'])+'\n')
                for domain in entry.topology:
                    f.write(','.join(map(str,[domain[1],domain[2],domain[0]]))+'\n')