# Notes de format pour Stat/Wiley

Les consignes de *Stat* indiquent que les articles doivent rester compacts, avec une longueur maximale d'environ 10 pages de texte, hors références, figures et tableaux. Source vérifiée le 2026-06-11 : <https://onlinelibrary.wiley.com/page/journal/20491573/homepage/forauthors.html>.

Wiley recommande son gabarit LaTeX lorsque pertinent et précise qu'une soumission LaTeX doit fournir un PDF compilé avec le fichier `.tex` principal et les fichiers de support, notamment figures et bibliographie. Source vérifiée le 2026-06-11 : <https://authors.wiley.com/author-resources/Journal-Authors/Prepare/latex-template.html>.

Ce projet n'inclut pas le paquet officiel `WileyDesign.zip`. Le manuscrit est donc préparé avec une classe `article` standard, une bibliographie BibTeX et un minimum de macros. Cette version est adaptée à une première soumission si le système accepte un PDF avec sources, et pourra être transposée dans le gabarit Wiley si l'interface l'exige.

Paquet local préparé :

1. `main.pdf`;
2. `main.tex`;
3. `references.bib`;
4. `figures/*.pdf` et `figures/*.png`;
5. `tables/simulation_table.tex`;
6. `supporting_information.pdf` et `supporting_information.tex`;
7. `numerics/` avec code, tests, sorties CSV/JSON et instructions.

Avant dépôt final :

1. Confirmer dans le système de soumission que le PDF avec sources LaTeX standard est accepté pour la première soumission.
2. Si le système demande le gabarit Wiley, copier le contenu de `main.tex` dans le gabarit sans modifier les énoncés mathématiques.
3. Vérifier le nombre de pages après toute conversion au gabarit officiel.
4. Confirmer les déclarations auteur : originalité, non-soumission ailleurs, financement, conflits d'intérêt, disponibilité du code.
