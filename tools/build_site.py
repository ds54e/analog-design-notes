"""MIT: copy an explicit static resource set and record the served revision.

No simulator, model checkout, scientific package or acquisition is imported.
"""
import hashlib,json,os,re,shutil,zipfile
from pathlib import Path

root=Path(__file__).resolve().parents[1]
release=json.loads((root/'RELEASE.json').read_text());slug=release['study'];example=root/'examples'/slug;site=root/'site';web=site/'studies'/slug
for name in release['example_files']:
    p=Path(name);assert not p.is_absolute() and '..' not in p.parts
    assert (example/p).is_file(),name
    assert hashlib.sha256((example/p).read_bytes()).hexdigest()==release['example_sha256'][name],name
manifest=json.loads((example/'data/manifest.json').read_text())
for name,h in manifest['files'].items():assert hashlib.sha256((example/'data'/name).read_bytes()).hexdigest()==h,name
for name in release['web_resources']:
    target=web/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(example/name,target)
licenses=site/'licenses';licenses.mkdir(exist_ok=True)
license_files=['LICENSE.md','THIRD_PARTY_NOTICES.md']+['LICENSES/'+p.name for p in sorted((root/'LICENSES').iterdir()) if p.is_file()]
for name in license_files:
    target=licenses/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(root/name,target)
downloads=site/'downloads';downloads.mkdir(exist_ok=True)
with zipfile.ZipFile(downloads/'gain-study-v1.0.zip','w',compression=zipfile.ZIP_DEFLATED) as z:
    for name in release['example_files']:z.write(example/name,arcname='example/'+name)
    for name in license_files:z.write(root/name,arcname='example/licenses/'+name)
revision=os.environ.get('GITHUB_SHA',os.environ.get('ADN_PREVIEW_REVISION','CANDIDATE'))
assert revision=='CANDIDATE' or re.fullmatch('[0-9a-f]{40}',revision)
(site/'revision.json').write_text(json.dumps(dict(commit=revision,study=slug,version=release['version']),indent=2)+'\n')
(site/'revision-header.html').write_text('<meta name="adn-source-commit" content="'+revision+'">\n')
(site/'revision-footer.html').write_text('<div id="adn-build">Study '+release['version']+' · source revision '+revision+'</div>\n')
print('Prepared only selected static resources, example download and revision metadata.')
