"""Fetch original inputs and verify their hashes; run from package root."""
from pathlib import Path
import gzip, hashlib, json, shutil, urllib.request
root = Path('genomics')
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 return h.hexdigest()
for item in json.loads((root/'source_manifest.json').read_text()):
 target=root/item['destination']; target.parent.mkdir(parents=True,exist_ok=True)
 if target.exists():
  if sha(target)!=item['sha256']:raise ValueError('Existing file differs: '+str(target))
  print('Verified',target);continue
 compressed=item.get('compression')=='gzip'
 download=target.with_name(target.name+'.download')
 request=urllib.request.Request(item['url'],headers={'User-Agent':'R (4.5.0 aarch64-apple-darwin20)'})
 with urllib.request.urlopen(request,timeout=300) as response,download.open('wb') as f:
  shutil.copyfileobj(response,f)
 if compressed:
  expanded=target.with_name(target.name+'.expanded')
  with gzip.open(download,'rb') as source,expanded.open('wb') as dest:shutil.copyfileobj(source,dest)
  download.unlink();download=expanded
 if sha(download)!=item['sha256']:raise ValueError('Downloaded file differs: '+str(target))
 download.rename(target);print('Downloaded and verified',target)
