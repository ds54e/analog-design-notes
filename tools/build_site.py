"""MIT: copy an explicit static resource set and record the served revision.

No simulator, model checkout, scientific package or acquisition is imported.
"""
import hashlib,json,os,re,shutil
from pathlib import Path

root=Path(__file__).resolve().parents[1]
release=json.loads((root/'RELEASE.json').read_text());site=root/'site'
downloads=site/'downloads';downloads.mkdir(exist_ok=True)
studies=[release]+release.get('additional_studies',[])
packages=studies+release.get('additional_packages',[])
for study in packages:
    slug=study['study'];example=root/'examples'/slug;web=site/'studies'/slug
    for name in study['example_files']:
        p=Path(name);assert not p.is_absolute() and '..' not in p.parts
        assert (example/p).is_file(),name
        assert hashlib.sha256((example/p).read_bytes()).hexdigest()==study['example_sha256'][name],name
    manifest=json.loads((example/'data/manifest.json').read_text())
    for name,h in manifest['files'].items():assert hashlib.sha256((example/'data'/name).read_bytes()).hexdigest()==h,name
    for name in study['web_resources']:
        target=web/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(example/name,target)
    archive=root/study['archive_file']
    assert hashlib.sha256(archive.read_bytes()).hexdigest()==study['archive_sha256']
    shutil.copyfile(archive,downloads/study['download'])
for name,expected in release.get('learning_calculations',{}).get('resource_sha256',{}).items():
    path=Path(name);assert not path.is_absolute() and '..' not in path.parts
    assert hashlib.sha256((site/path).read_bytes()).hexdigest()==expected,name
licenses=site/'licenses';licenses.mkdir(exist_ok=True)
license_files=['LICENSE.md','THIRD_PARTY_NOTICES.md']+['LICENSES/'+p.name for p in sorted((root/'LICENSES').iterdir()) if p.is_file()]
for name in license_files:
    target=licenses/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(root/name,target)
revision=os.environ.get('GITHUB_SHA',os.environ.get('ADN_PREVIEW_REVISION','CANDIDATE'))
assert revision=='CANDIDATE' or re.fullmatch('[0-9a-f]{40}',revision)
(site/'revision.json').write_text(json.dumps(dict(commit=revision,study=release['study'],version=release['version'],
    collection_version=release.get('collection_version'),learning_edition=release.get('learning_edition'),
    studies=[dict(study=s['study'],version=s['version'],research_tag=s['research_tag']) for s in studies],
    additional_packages=[dict(package=s['study'],version=s['version'],tag=s['research_tag']) for s in release.get('additional_packages',[])]),indent=2)+'\n')
(site/'revision-header.html').write_text('<meta name="adn-source-commit" content="'+revision+'">\n')
(site/'revision-footer.html').write_text('<div id="adn-build">Collection '+release.get('collection_version',release['version'])+' · source revision '+revision+'</div>\n')
print('Prepared only selected static resources, example download and revision metadata.')
