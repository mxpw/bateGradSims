# bateGradSims

Assessing sampling strategies effects on Bateman gradients estimates

## Install

devtools is needed

```bash
install.package('devtools')
devtools::install_github("mxpw/bateGradSims")
```

## As usual


Functions are more or less documented. Use classic R style to see doc.

```bash
library(bateGradSims)
?ms_obs
```

## Changes logs

### TODO

- [ ] Proposition pour lier nombre de gamètes et compétitivité pour les mâles (covariance positive, négative, nulle)
  - [ ] Utiliser distributions conditionnelles ? (i.e. ~P(Y|X))
- [ ] "Corrélation de parentalité" i.e. probabilité que deux graines tirées aléatoirement soient du même père (voir Dorken & Perry, 2017) 
  - [ ] Version de base
  - [ ] Version avec correction pour tailles de pop. 
- [ ] Implémenter la gestion auto. des plots des gradients
	- [ ] Ajouter annotations (e.g., valeur de $\beta$, éventuellement $\Delta(\beta_{true}-\beta_{est})$)
	- [ ] Attention - en l'état, ce ne sont pas les $\beta$ estimés qui sont représentés
- [ ] Ajout statistiques descriptives sur MS/RS (moyennes et variances) - sur les échantillons ainsi que sur les RS/MS vrais
- [x] Bootstrap des tirages
	- [ ] Revoir la fonction sampling pour intégrer par défaut le groundtruth() (et ne pas le répliquer...)
- [x] Paternity share 
- [ ] Extraire la partie qui créer la matrice des alpha de dirichlet de la fonction pollen export - sera plus flexible (permet de définir différentes méthodes pour créer la dite matrice voir d'en spécifier une à la main)
  - [ ] TODO - définir method
- [ ] Voir pour factoriser les bouts de codes communs dans les fonctions sampling_XXX()
- [ ] Notes explicatives
  - [ ] Dirichlet & Multinomiale (pollen_export)
  - [ ] autres trucs


### Checks
- [ ] Cas où pas d'ovule fécondé pour une femelle (NULL/integer(0)/NAN/0) ? (problèmatique gestion du code, et sens biologique)
- [ ] Question des arrondis - notamment pro-rata ; 10% sur 4 ovules -> 0.4 => 0 ou 1 ?
- [ ] Check valeurs par défaut pour ms/rs (NA vs. 0)
- [ ] Sampling chelou (x/mean(x)) - voir fonction sampling() + commentaire à la fin
- [ ] Check des group_by avec replicats 
  - [ ] pour le sampling
  - [ ] pour le fit gradients
- [ ] Revoir toute la doc >_<