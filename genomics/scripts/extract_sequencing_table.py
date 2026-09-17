"""Extract numeric supplementary tables, retaining the original genome build."""
from pathlib import Path
import csv,re
from pypdf import PdfReader
root=Path('genomics/data/genomics')
rows=[];cell=None
pattern=re.compile(r'^\s*(\d+|X|Y)\s+([\d,]+)\s+([\d,]+)\s+([\d.]+)\s+([\d.Ee+\-]+)\s+([\d,]+)\s+([\d,]+)\s*$')
for page in PdfReader(root/'Chiang2009_sequencing_segments_hg18.pdf').pages:
 text=page.extract_text()
 for label in ['HCC1954','HCC1143','NCI-H2347']:
  if ('Segmentation of '+label) in text:cell=label
 for line in text.splitlines():
  match=pattern.match(line)
  if not match:continue
  assert cell is not None
  chrom,start,end,ratio,pval,normal,tumor=match.groups()
  rows.append([cell,chrom,int(start.replace(',','')),int(end.replace(',','')),float(ratio),float(pval),int(normal.replace(',','')),int(tumor.replace(',',''))])
assert sum(row[0]=='HCC1143' for row in rows)==407
with (root/'Chiang2009_sequencing_segments_hg18.csv').open('w',newline='') as f:
 writer=csv.writer(f);writer.writerow(['cell_line','chromosome','start','end','copy_ratio','p_value','normal_count','tumor_count']);writer.writerows(rows)
print('Extracted',len(rows),'segments from the publisher tables (hg18).')
