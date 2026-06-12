# Projet Stat/Wiley - Corrective hidden states in misspecified HMMs

Ce depot contient le manuscrit, le supplement, le code numerique reproductible
et les fichiers de preparation pour une soumission a *Stat* de Wiley.

## Contenu

```text
article/
  main.tex                  Manuscrit LaTeX compact
  main.pdf                  PDF compile du manuscrit
  references.bib            Bibliographie BibTeX
  figures/                  Figures numeriques
  tables/                   Tables LaTeX generees
supplement/
  supporting_information.tex
  supporting_information.pdf
  tables/                   Tables du supplement
numerics/
  scripts/run_review_simulation.py
  src/                      Simulation, EM HMM, diagnostics et figures
  tests/                    Tests numeriques
  output/                   Resultats CSV/JSON/LaTeX generes
tools/
  validate_submission.py
  build_submission_package.py
submission/
  stat_hmm_submission_package.zip
```

## Manuscrit

Le manuscrit est en anglais. Il est volontairement compact et centre sur le
Theoreme 5.1: sous des conditions suffisantes de separation, un HMM raffine
reproduit exactement la loi observee, tandis qu'un HMM non raffine avec le
nombre d'etats agreges reste separe en taux de Kullback-Leibler.

Le supplement contient les details de l'EM a emissions melange, les resultats
des 200 replications par taille d'echantillon et les diagnostics graphiques
additionnels.

## Reproduction numerique

Depuis la racine du projet:

```bash
python3 numerics/scripts/run_review_simulation.py
python3 -m pytest numerics/tests -q
```

Le script genere:

- les figures PDF/PNG dans `article/figures/`;
- la table principale dans `article/tables/simulation_table.tex`;
- la table du supplement dans `supplement/tables/review_summary_table.tex`;
- les resultats CSV/JSON dans `numerics/output/`.

Aucune donnee externe n'est utilisee.

## Compilation

Depuis `article/`:

```bash
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

Depuis `supplement/`:

```bash
pdflatex supporting_information.tex
pdflatex supporting_information.tex
```

## Validation locale

Depuis la racine du projet:

```bash
python3 tools/validate_submission.py
python3 tools/build_submission_package.py
```

Le validateur verifie les fichiers attendus, les marqueurs non publics, les
sorties de 200 replications par taille d'echantillon et un ajustement numerique
court incluant le HMM a emissions melange.

## Paquet de soumission

`tools/build_submission_package.py` reconstruit
`submission/stat_hmm_submission_package/` et
`submission/stat_hmm_submission_package.zip` avec le manuscrit, le supplement,
les sources LaTeX, les figures, les tables, le code numerique et les sorties
reproductibles.

Avant depot final, confirmer dans le systeme de soumission les declarations
auteur: originalite, non-soumission ailleurs, financement, conflits d'interet
et disponibilite du code.
