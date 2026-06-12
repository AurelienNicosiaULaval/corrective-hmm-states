# Stat/Wiley format notes

The *Stat* author guidelines indicate that articles should remain compact, with
an approximate maximum length of 10 text pages, excluding references, figures
and tables. Source checked on 2026-06-11:
<https://onlinelibrary.wiley.com/page/journal/20491573/homepage/forauthors.html>.

Wiley recommends its LaTeX template when appropriate and states that a LaTeX
submission should include a compiled PDF together with the main `.tex` file and
supporting files, including figures and bibliography. Source checked on
2026-06-11:
<https://authors.wiley.com/author-resources/Journal-Authors/Prepare/latex-template.html>.

This project does not include the official `WileyDesign.zip` template package.
The manuscript is therefore prepared with a standard `article` class, a BibTeX
bibliography, and a minimal macro set. This version is suitable for an initial
submission if the system accepts a PDF with source files, and it can be moved
into the Wiley template if the submission interface requires it.

Prepared local package:

1. `main.pdf`;
2. `main.tex`;
3. `references.bib`;
4. `figures/*.pdf` and `figures/*.png`;
5. `tables/simulation_table.tex`;
6. `supporting_information.pdf` and `supporting_information.tex`;
7. `numerics/` with code, tests, CSV/JSON outputs and instructions.

Before final upload:

1. Confirm in the submission system that a standard LaTeX source package is
   accepted for first submission.
2. If the system requires the Wiley template, copy the content of `main.tex`
   into the template without changing the mathematical statements.
3. Check the page count after any conversion to the official template.
4. Confirm author declarations: originality, no concurrent submission
   elsewhere, funding, conflicts of interest and code availability.
