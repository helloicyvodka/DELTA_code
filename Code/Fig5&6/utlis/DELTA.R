# required options:
# -treeS: TreeS file path
# -treeT: TreeT file path
# -cost: Cost file path. the cost file contains the score for different types of leaves


# optional parameters:
# -method: l or g. g for global alignment; l for local alignment. default: g.
# -max_target: target num for l, local alignment. default: 1
# -test: testNum to calculate p-value. If testNum <=2, do not output p-value. default 0
# -outfile: output file path. default: TreeS file path + l or g, based on -method
# -all: T or F. output as much information as possible. default F;
# -prune: pruneScore is the punish for pruning one leaf. default 1

DELTA.address <- "/mnt/data/home/phil/acting/treeComparison/code/DELTA/bin/Release/DELTA"



DELTA <- function(treeS,
                  treeT,
                  cost,
                  outfile,
                  method="g",
                  max_target=1,
                  test=0,
                  all="F",
                  prune=1,
                  DELTA.address="/mnt/data/home/phil/acting/treeComparison/code/DELTA/bin/Release/DELTA"){
  
  system(paste0(DELTA.address,
                " ",
                "-treeS",
                " ",
                treeS,
                " ",
                "-treeT",
                " ",
                treeT,
                " ",
                "-cost",
                " ",
                cost,
                " ",
                "-max_target",
                " ",
                max_target,
                " ",
                "-method",
                " ",
                method,
                " ",
                "-outfile",
                " ",
                outfile,
                " ",
                "-test",
                " ",
                test,
                " ",
                "-all",
                " ",
                all,
                " ",
                "-prune",
                " ",
                prune
                ),
  intern = F)
  
}
