import IMP
import IMP.pmi
import macros_e29
import sys

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

print inputs

is_mpi = inputs.is_mpi                                          # mpi enabled
load_distance_matrix_file = inputs.load_distance_matrix_file    # skip the matrix calcuklation and read the precalculated matrix
nbestscoringmodels = int(inputs.nbestscoringmodels)             # number of models to be clustered
nclusters = int(inputs.nclusters)                               # number of clusters to be used by kmeans algorithm
prefiltervalue = float(inputs.prefiltervalue)                   # prefilter the models by score


#####################################################
# initialize the macro
#####################################################
model=IMP.Model()

mergedirectories= []
mergedirectories.append("./all_models.999")


mc=macros_e29.AnalysisReplicaExchange0(model,
                                        stat_file_name_suffix="stat",
                                        merge_directories=mergedirectories,
                                        global_output_directory="./")

feature_list=["Total_Score",
              "ConnectivityRestraint_None",
              "CrossLinkingMassSpectrometryRestraint_Data_Score",
              "ExcludedVolumeSphere_all",
              "ExcludedVolumeSphere_ecm29"
              ]

# Dictionary of densities to be calculated
# the key is the name of the file and the value if the selection
# example:
#              {"med17-CTD":[(200,300,"med17")],"med17-CTD.med14":[(200,300,"med17"),"med14"]   }

reduced_density_dict={"Rpt1":["Rpt1"], "Rpt2":["Rpt2"], "Rpt3":["Rpt3"], 
                      "Rpt4":["Rpt4"], "Rpt5":["Rpt5"], "Rpt6":["Rpt6"],
                      "RptS":["Rpt1", "Rpt2", "Rpt3", "Rpt4", "Rpt5", "Rpt6"],
                      "Rpn1":["Rpn1"], "Rpn2":["Rpn2"], "Rpn3":["Rpn3"], "Rpn5":["Rpn5"],
                      "Rpn6":["Rpn6"], "Rpn7":["Rpn7"], "Rpn8":["Rpn8"], "Rpn9":["Rpn9"],
                      "Rpn10":["Rpn10"], "Rpn11":["Rpn11"], "Rpn12":["Rpn12"], "Rpn15":["Rpn15"],
                      "RpnS":["Rpn1", "Rpn2", "Rpn3", "Rpn5", "Rpn6", "Rpn7", "Rpn8", "Rpn9", "Rpn10", "Rpn11", "Rpn12", "Rpn15"],
                      
                      "ecm29_352":[(1,504,"ecm29")], "ecm29_689":[(505, 1845, "ecm29")], 
                      "ecm29S":[(1,504,"ecm29"), (505, 1845, "ecm29")],
                      "Whole":["Rpt1", "Rpt2", "Rpt3", "Rpt4", "Rpt5", "Rpt6", "Rpn1", "Rpn2", "Rpn3", "Rpn5", "Rpn6", "Rpn7", "Rpn8", "Rpn9", "Rpn10", "Rpn11", "Rpn12", "Rpn15", (1,504,"ecm29"),(505, 1845, "ecm29")]
                      }                
# list of component names needed to calculate the RMSD for the clustering
components_names={"ecm29_352":(1,504,"ecm29"), 
                  "ecm29_689":(505, 1845, "ecm29")}

mc.clustering("Score_None",  # don't change, field where to find the score
              "rmf_file",                          # don't change, field where to find the path for the rmf_file
              "rmf_frame_index",                   # don't change, field for the frame index
              number_of_best_scoring_models=nbestscoringmodels,   # number of models to be clustered
              alignment_components=None,         # don't change, (list of proteins you want to use for structural alignment
              rmsd_calculation_components=components_names,  # list of proteins used to calculated the rmsd
              distance_matrix_file="distance.rawmatrix.pkl", # save the distance matrix
              outputdir="kmeans_"+str(nbestscoringmodels)+"_"+str(nclusters)+"/",  # directory name for the clustering
              feature_keys=feature_list,                     # extract these fields from the stat file
              load_distance_matrix_file=load_distance_matrix_file,                # skip the matrix calcuklation and read the precalculated matrix
              skip_clustering=False,                         # skip clustering
              display_plot=True,                            # display the heat map plot of the distance matrix
              exit_after_display=False,                      # exit after having displayed the distance matrix plot
              get_every=1,                                   # skip structures for faster computation
              number_of_clusters=nclusters,                  # number of clusters to be used by kmeans algorithm
              voxel_size=5.0,                                # voxel size of the mrc files
              density_custom_ranges=reduced_density_dict)    # setup the list of densities to be calculated

