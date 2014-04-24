import argparse

parser = argparse.ArgumentParser() 
parser.add_argument("--launch-onto-cluster",help="launches the whole thing on the cluster")
args = parser.parse_args()  
print args.launch_onto_cluster
