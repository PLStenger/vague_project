# vague_project
Projet VAGUE - nouVeaux usages et Adaptations du microbiome marin face aux chanGements climatiques et environnementaUx observEs en zone Pacifique Tropicale - co-portés par Véronique Anton (CNRS/IRD) et Delphine Dissard-Field (IRD)

### Installing pipeline :

First, open your terminal. Then, run these two command lines :

    cd -place_in_your_local_computer
    git clone https://github.com/PLStenger/vague_project.git

### Update the pipeline in local by :

    git pull

### Running the script with :

    time nohup bash 01_qiime2_preprocess.sh &> 01_qiime2_preprocess.out
    time nohup bash 01_qiime2_preprocess_PA_KNS.sh &> 01_qiime2_preprocess_PA_KNS.out
    time nohup bash 02_diversity_analysis.sh &> 02_diversity_analysis.out
