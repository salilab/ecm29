import ihm
import os

# If we're running from an SGE job, override the from_pubmed_id() function
# to return a cached value, since we don't have network access (needed to
# query PubMed directly)

def mock_from_pubmed(cls, pubmed_id):
    return ihm.Citation(
            pmid=28821611,
            title='The proteasome-interacting Ecm29 protein disassembles the '
                  '26S proteasome in response to oxidative stress.',
            journal='J Biol Chem', volume=292, page_range=(16310,16320),
            year=2017, doi='10.1074/jbc.M117.803619', authors=[
                  'Wang X', 'Chemmama IE', 'Yu C', 'Huszagh A', 'Xu Y',
                  'Viner R', 'Block SA', 'Cimermancic P', 'Rychnovsky SD',
                  'Ye Y', 'Sali A', 'Huang L'])

if 'JOB_ID' in os.environ:
    ihm.Citation.from_pubmed_id = classmethod(mock_from_pubmed)
