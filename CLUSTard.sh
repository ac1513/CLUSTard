#!/bin/bash

while getopts ":w:" opt; do
  case $opt in
    w)
      wf=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

if ! command -v snakemake &> /dev/null
then
    echo " ";
    echo " ===========================================";
    echo "   ____ _    _   _ ____ _____             _ ";
    echo "  / ___| |  | | | / ___|_   _|_ _ _ __ __| |";
    echo " | |   | |  | | | \___ \ | |/ _\` | '__/ _\` |";
    echo " | |___| |__| |_| |___) || | (_| | | | (_| |";
    echo "  \____|_____\___/|____/ |_|\__,_|_|  \__,_|";
    echo "                                            ";
    echo " ===========================================";
    echo "\nSnakemake could not be found, Please ensure Snakemake (>v5.3) is installed before using the pipeline.\n"
    exit
fi

if [[ -n $wf ]] && [[ $wf == "both" ]];
then
  echo " ";
  echo " ===========================================";
  echo "   ____ _    _   _ ____ _____             _ ";
  echo "  / ___| |  | | | / ___|_   _|_ _ _ __ __| |";
  echo " | |   | |  | | | \___ \ | |/ _\` | '__/ _\` |";
  echo " | |___| |__| |_| |___) || | (_| | | | (_| |";
  echo "  \____|_____\___/|____/ |_|\__,_|_|  \__,_|";
  echo "                                            ";
  echo " ===========================================";
  echo "  Welcome to CLUSTard, binning and analysis ";
  echo "  workflows will be run as you specified   ";
  echo "  both! 3... 2... 1...";
  echo " ";
  { snakemake -n --snakefile workflow/Snakefile_bin &&
  echo " Thank you for using CLUSTard, your metagenome bins are in the directory: " &&
  echo "         results/clusters/                  " &&
  echo " Have a good day!"
} || echo "\nError: Something went wrong, please see above for error message. \n "
  { snakemake -n --snakefile workflow/Snakefile_analysis &&
  echo " Thank you for using CLUSTard, your metagenome bins are in the directory: " &&
  echo "         results/clusters/                  " &&
  echo " Have a good day!"
} || echo "\nError: Something went wrong, please see above for error message. \n"

elif [[ -n $wf ]] && [[ $wf == "bin" ]] || [[ $wf == "binning" ]];
then
  echo " ";
  echo " ===========================================";
  echo "   ____ _    _   _ ____ _____             _ ";
  echo "  / ___| |  | | | / ___|_   _|_ _ _ __ __| |";
  echo " | |   | |  | | | \___ \ | |/ _\` | '__/ _\` |";
  echo " | |___| |__| |_| |___) || | (_| | | | (_| |";
  echo "  \____|_____\___/|____/ |_|\__,_|_|  \__,_|";
  echo "                                            ";
  echo " ===========================================";
  echo "  Welcome to CLUSTard, the binning workflow ";
  echo "  will be run as you specified $wf          ";
  echo "  3... 2... 1...";
  echo " ";
  { snakemake -n --snakefile workflow/Snakefile_bin &&
  echo " Thank you for using CLUSTard, your metagenome bins are in the directory: " &&
  echo "         results/clusters/                  " &&
  echo " Have a good day!"
} || echo "\nError: Something went wrong, please see above for error message. \n"
elif [[ -n $wf ]] && [[ $wf == "analysis" ]] || [[ $wf == "ana" ]];
then
  echo " ";
  echo " ===========================================";
  echo "   ____ _    _   _ ____ _____             _ ";
  echo "  / ___| |  | | | / ___|_   _|_ _ _ __ __| |";
  echo " | |   | |  | | | \___ \ | |/ _\` | '__/ _\` |";
  echo " | |___| |__| |_| |___) || | (_| | | | (_| |";
  echo "  \____|_____\___/|____/ |_|\__,_|_|  \__,_|";
  echo "                                            ";
  echo " ===========================================";
  echo "  Welcome to CLUSTard, the analysis workflow ";
  echo "  will be run as you specified $wf          ";
  echo "  The location of your metagenome bins is:  ";
  echo "  TBD";
  echo "  3... 2... 1...";
  echo " ";
  { snakemake -n --snakefile workflow/Snakefile_bin &&
  echo " Thank you for using CLUSTard, your metagenome bins are in the directory: " &&
  echo "         results/clusters/                  " &&
  echo " Have a good day!"
} || echo "\nError: Something went wrong, please see above for error message. \n "

else [[ -n $wf ]];
  echo " $wf is not a valid option, please specify workflow from";
  echo " 'both', 'binning' or 'analysis'  ";
fi
