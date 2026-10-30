# Projet FEM - Analyse de Plasticité

Ce projet implémente un solveur d'analyse par éléments finis (FEA) incrémental-itératif pour la plasticité de von Mises. Il modélise à la fois l'écrouissage isotrope (expansion de la surface de charge) et l'écrouissage cinématique (translation de la surface de charge via la contrainte inverse).

## Structure des Fichiers

### Code Source Fortran
* **`main.f90`** : Programme principal orchestrant la boucle incrémentale-itérative globale et les itérations de Newton.
* **`types.f90`, `fedata.f90`, `fea.f90`** : Structures de données, paramètres globaux et routines FEA.
* **`plane42.f90`** : Élément fini 2D à 4 nœuds et assemblage de la matrice de rigidité tangente.
* **`link1.f90`** : Éléments de type barre unidimensionnelle.
* **`numeth.f90`, `processor.f90`** : Méthodes numériques et traitement des étapes.

### Fichiers d'Entrée et Données
* **`test.txt`** : Fichier d'entrée avec la géométrie et les incréments de charge.
* **`results.dat`** : Fichier de sortie principal (déplacements, forces).
* **`stress_log.dat`** : Journal des états de contraintes.
* **`gauss_stress_history.dat`** : Historique des contraintes aux points de Gauss.

### Scripts de Post-Traitement (MATLAB)
* **`plot_stress_strain.m`** : Visualisation des courbes contrainte-déformation.
* **`plasticity_scale.m`** : Distribution spatiale de la déformation plastique.

## Utilisation
1. Compiler le projet :
   ```bash
   gfortran types.f90 fedata.f90 numeth.f90 plane42.f90 link1.f90 fea.f90 processor.f90 main.f90 -o fem_solver