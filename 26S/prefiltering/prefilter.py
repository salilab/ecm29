import IMP
import IMP.pmi
#import IMP.pmi.macros
import sys,os
import glob
#####################################################
# Parsing parameter inputs
#####################################################
import argparse

parser = argparse.ArgumentParser(description='generate clusters of the RMF solutions')
parser.add_argument('-mpi', action="store", dest="is_mpi", help="mpi enabled")
parser.add_argument('-preload', action="store", dest="load_distance_matrix_file", help="skip the matrix calcuklation and read the precalculated matrix")
parser.add_argument('-nmods', action="store", dest="nbestscoringmodels", help="number of models to be clustered")
parser.add_argument('-nclusters', action="store", dest="nclusters", help="number of clusters to be used by kmeans algorithm")
parser.add_argument('-prefilter', action="store", dest="prefiltervalue", help="prefilter the models by score")
parser.add_argument('-dir', action="store", dest="is_dir", help="input directory")
inputs = parser.parse_args()

# mpi enabled
if (inputs.is_mpi=="True") or (inputs.is_mpi=="true") or (inputs.is_mpi=="Yes") or (inputs.is_mpi=="yes") :
    inputs.is_mpi = True
else:
    inputs.is_mpi = False

# skip the matrix calcuklation and read the precalculated matrix
if (inputs.load_distance_matrix_file=="True") or (inputs.load_distance_matrix_file=="true") or (inputs.load_distance_matrix_file=="Yes") or (inputs.load_distance_matrix_file=="yes") :
    inputs.load_distance_matrix_file = True
else:
    inputs.load_distance_matrix_file = False

# number of models to be clustered
if inputs.nbestscoringmodels==None:
    inputs.nbestscoringmodels = 500

# number of clusters to be used by kmeans algorithm
if inputs.nclusters==None:
    inputs.nclusters = 1

# prefilter the models by score
if inputs.prefiltervalue==None:
    inputs.prefiltervalue = 565.0
    
if inputs.is_dir==None:
    exit
print inputs

is_mpi = inputs.is_mpi                                          # mpi enabled
load_distance_matrix_file = inputs.load_distance_matrix_file    # skip the matrix calcuklation and read the precalculated matrix
nbestscoringmodels = int(inputs.nbestscoringmodels)             # number of models to be clustered
nclusters = int(inputs.nclusters)                               # number of clusters to be used by kmeans algorithm
prefiltervalue = float(inputs.prefiltervalue)                   # prefilter the models by score
is_dir = inputs.is_dir

#####################################################
# initialize the macro
#####################################################
import macros_e29
model=IMP.Model()


#for dir in os.listdir('/salilab/park3/peterc/Lan_XL/ecm29/pmi_modeling_try3/outputs/output_?'):
#    dir=os.path.join("/salilab/park3/peterc/Lan_XL/ecm29/pmi_modeling_try3/outputs",dir)

print(is_dir)
mergedirectories= []
mergedirectories.append(is_dir)    
mc=macros_e29.AnalysisReplicaExchange0(model,
                                       stat_file_name_suffix="stat",     # don't change
                                       merge_directories=mergedirectories,
                                       global_output_directory='output')


feature_list=["Total_Score",
              "ConnectivityRestraint_None",
              "CrossLinkingMassSpectrometryRestraint_Data_Score",
              "ExcludedVolumeSphere_all",
              "ExcludedVolumeSphere_ecm29"
              ]

reduced_density_dict=None

components_names={"ecm29_352":(1,504,"ecm29"), 
                  "ecm29_689":(505, 1845, "ecm29")}

mc.clustering("Total_Score",  
              "rmf_file",                          
              "rmf_frame_index",                  
              prefiltervalue=prefiltervalue,          
              number_of_best_scoring_models=nbestscoringmodels,
              alignment_components=components_names,
              rmsd_calculation_components=components_names,
              distance_matrix_file="distance.rawmatrix.pkl",
              outputdir="kmeans_"+str(nbestscoringmodels)+"_"+is_dir.split("/")[-1]+"/", 
              feature_keys=feature_list,
              load_distance_matrix_file=load_distance_matrix_file,
              skip_clustering=True,
              display_plot=True,
              exit_after_display=False,
              get_every=1,
              number_of_clusters=nclusters,
              voxel_size=5.0,
              density_custom_ranges=reduced_density_dict)
