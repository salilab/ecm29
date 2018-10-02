from __future__ import print_function
import IMP
import RMF
import IMP.atom
import IMP.core
import IMP.algebra
import IMP.container
import IMP.rmf
import os, sys
import ihm
import IMP.pmi
import IMP.pmi.mmcif
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking


#---------------------------
# Define Input Files
#---------------------------
datadirectory = "../data/"
topology_file = datadirectory+"topology.txt"

#--------------------------
# Set MC Sampling Parameters
#--------------------------
num_frames = 7500
if '--test' in sys.argv: num_frames=20
num_mc_steps = 10

#--------------------------
# Create movers
#--------------------------

# rigid body movement params
rb_max_trans = 1.00
rb_max_rot = 0.01
# flexible bead movement
bead_max_trans = 2.00


#--------------------------------
# Build the Model Representation
#--------------------------------

# Initialize model
m = IMP.Model()

# Create list of components from topology file
topology = IMP.pmi.topology.TopologyReader(topology_file,
                                           pdb_dir = datadirectory,
                                           fasta_dir = datadirectory,
                                           )
domains = topology.get_components()
print('#'*10,domains)

bs = IMP.pmi.macros.BuildSystem(m)
bs.dry_run = '--dry-run' in sys.argv

if '--mmcif' in sys.argv:
    # Record the modeling protocol to an mmCIF file
    po = IMP.pmi.mmcif.ProtocolOutput(open('ecm29.cif', 'w'))
    bs.system.add_protocol_output(po)
    po.system.title = ('The proteasome-interacting Ecm29 protein disassembles '
                       'the 26S proteasome in response to oxidative stress')
    # Add publication
    po.system.citations.append(ihm.Citation.from_pubmed_id(28821611))

bs.add_state(topology)
representation, dof = bs.execute_macro(max_rb_trans=rb_max_trans, 
                                       max_rb_rot=rb_max_rot, 
                                       max_bead_trans=bead_max_trans)
#representation.shuffle_configuration(50)

#--------------------------
# Define Degrees of Freedom
#--------------------------

# Add default mover parameters to simulation
outputobjects = [] # reporter objects (for stat files)
sampleobjects = [] # sampling objects

#-----------------------------------
# Define Scoring Function Components
#-----------------------------------

# Here we are defining a number of restraints on our system.
#  For all of them we call add_to_model() so they are incorporated into scoring
#  We also add them to the outputobjects list, so they are reported in stat files

# Excluded Volume Restraint
#  To speed up this expensive restraint, we operate it at resolution 20

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print("ilan0", sf.evaluate(False))

ecm29 = []
p26ps = []
crs = []
moldict = bs.get_molecules()[0]
for molname in moldict:
   for mol in moldict[molname]:      
      cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol, scale=2.0)
      cr.add_to_model()
      outputobjects.append(cr)
      sampleobjects.append(cr)
      crs.append(cr)

      if 'ecm29' in molname:
         atomic = mol.get_atomic_residues()
         dof.create_rigid_body(mol[0:503] & atomic,
                                 max_trans=rb_max_trans,
                                 max_rot=rb_max_rot,
                                 resolution='all')
         dof.create_rigid_body(mol[504:1844] & atomic,
                               max_trans=rb_max_trans,
                               max_rot=rb_max_rot,
                               resolution='all')
         dof.create_flexible_beads(mol.get_non_atomic_residues(),
                                   max_trans=bead_max_trans,
                                   resolution=10)
         ecm29.append(mol)
      elif 'ecm29' not in molname:
         p26ps.append(mol)
         dof.create_flexible_beads(mol.get_non_atomic_residues(),
                                   max_trans=bead_max_trans,
                                   resolution=10)

print(ecm29)
print(p26ps)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print("ilan0", sf.evaluate(False))
ev1 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects = ecm29,
                                                              resolution=10)
ev1.add_to_model()
ev1.set_label('ecm29')
outputobjects.append(ev1)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print("ilan1", sf.evaluate(False))

ev2 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects = ecm29,
                                                              other_objects = p26ps,
                                                              resolution=10)

ev2.add_to_model()
ev2.rs.set_weight(1.0)
ev2.set_label('all')
outputobjects.append(ev2)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print("ilan2", sf.evaluate(False))

# Crosslinks - dataset 1
#  To use this restraint we have to first define the data format
#  Here assuming that it's a CSV file with column names that may need to change
#  Other options include the linker length and the slope (for nudging components together)
kw = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
kw.set_unique_id_key("id")
kw.set_protein1_key("prot1")
kw.set_protein2_key("prot2")
kw.set_residue1_key("res1")
kw.set_residue2_key("res2")
kw.set_id_score_key(None)
xldb = IMP.pmi.io.crosslink.CrossLinkDataBase(kw)
xldb.create_set_from_file(datadirectory+'xlinks.csv')

xl1 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    root_hier=representation,
    CrossLinkDataBase=xldb,
    length=21,
    label="Lan",
    resolution=1.0,
    slope=0.02)

xl1.rs.set_weight(25.0)
xl1.add_to_model()             # crosslink must be added to the model

sampleobjects.append(xl1) #crosslink restraint is storing a sampled particle
outputobjects.append(xl1)
dof.get_nuisances_from_restraint(xl1)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print("ilan3", sf.evaluate(False))

if not bs.dry_run:
    IMP.pmi.tools.shuffle_configuration(ecm29)
    dof.optimize_flexible_beads(100)

#--------------------------
# Monte-Carlo Sampling
#--------------------------
# This object defines all components to be sampled as well as the sampling protocol
mc1=IMP.pmi.macros.ReplicaExchange0(m,
                                    root_hier=representation,
                                    monte_carlo_sample_objects=dof.get_movers(),
                                    output_objects=outputobjects,
                                    crosslink_restraints=sampleobjects,
                                    monte_carlo_temperature=1.0,
                                    simulated_annealing=True,
                                    simulated_annealing_minimum_temperature=1.0,
                                    simulated_annealing_maximum_temperature=1.5,
                                    simulated_annealing_minimum_temperature_nframes=200,
                                    simulated_annealing_maximum_temperature_nframes=20,
                                    number_of_best_scoring_models=10,
                                    monte_carlo_steps=num_mc_steps,
                                    number_of_frames=num_frames,
                                    global_output_directory="output",
                                    test_mode=bs.dry_run)

# Start Sampling
mc1.execute_macro()

if '--mmcif' in sys.argv:
    # Correct number of output models to account for multiple runs
    protocol = po.system.orphan_protocols[-1]
    protocol.steps[-1].num_models_end = 3750000
    # Next, we filtered down to 109951 good scoring models
    analysis = ihm.analysis.Analysis()
    protocol.analyses.append(analysis)
    analysis.steps.append(ihm.analysis.FilterStep(
                            feature='energy/score',
                            num_models_begin=3750000, num_models_end=109951))

    # Finally, we found two clusters
    clusters = [{'rmf':'280_5415.rmf3'},
                {'rmf':'307_6895.rmf3'}]
    for ncluster, cluster in enumerate(clusters):
        with open('../Results/clustering/cluster.%d.all.txt' % ncluster) as fh:
            cluster['size'] = len(fh.readlines())

    analysis.steps.append(ihm.analysis.ClusterStep(
                            feature='RMSD', num_models_begin=109951,
                            num_models_end=sum(x['size'] for x in clusters)))

    for ncluster, cluster in enumerate(clusters):
        e = po._add_simple_ensemble(analysis.steps[-1],
                                    name="Cluster %d" % ncluster,
                                    num_models=cluster['size'],
                                    drmsd=None, num_models_deposited=1,
                                    localization_densities={},
                                    ensemble_file=None)
        # Add localization density for ecm29
        loc = ihm.location.OutputFileLocation(
               '../Results/localizations_densities/%d_ecm29.mrc' % (ncluster+1))
        den = ihm.model.LocalizationDensity(file=loc,
                                            asym_unit=po.asym_units['ecm29.0'])
        e.densities.append(den)

        # Add one output model - todo: doesn't currently work because RMF
        # topology doesn't match that here
#       rh = RMF.open_rmf_file_read_only('../Results/clustering/%s'
#                                        % cluster['rmf'])
#       IMP.rmf.link_hierarchies(rh, [representation])
#       IMP.rmf.load_frame(rh, RMF.FrameID(0))

#       model = po.add_mdoel(e.model_group)

    # Correct crosslinker type from "Lan" to DSSO
    for r in po.system.restraints:
        if hasattr(r, 'linker_type') and r.linker_type == 'Lan':
            r.linker_type = 'DSSO'

    po.flush()
